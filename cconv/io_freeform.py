# This is a Python version of the CASTEP io module - specifically the 'freeform' file reading
#
# It is heavily inspired in its interface to the original io module but the code is written from scratch and uses a more object-oriented approach

import re, collections, copy
from units import phys_units, default_units

class io_keyword_error(Exception):
    pass

class io_freeform_error(Exception):
    pass

class keyword(object):
    """Keyword class - defines a keyword to be used when parsing a freeform file
    
    Args:
    
       key(str): the actual keyword to be used in the file
       
       typ(str): written in the format X:Y, where X stands for the type of the variable associated with this keyword and Y for its 'level'. Possible values for X are:
       
       - S = String
       - I = Integer                                                       
       - R = Real                                                          
       - P = Physical                                                      
       - D = Defined                                                       
       - L = Boolean (Logical)                                             
       - V = Real Vector                                                   
       - W = Integer Vector                                                
       - B = Block
       
       Possible values for Y are:
       
       - B = Basic                                                         
       - I = Intermediate                                                  
       - E = Expert                                                        
       - D = Dummy (Unused)                                                
       
       description(str, optional): provides a helpful description of the meaning of the variable associated with the keyword
    
    Attrs:
        
        same as arguments
       
    """
    
    def __init__(self, key, typ, description=''):
        
        # Check that arguments are ok
        if type(key) is not str:
            raise io_keyword_error("Invalid key argument passed to keyword")
        else:
            self.key = key
        
        try:
            typmatch = re.findall('(S|I|R|P|D|L|V|W|B):(B|I|E|D)', typ.upper()) 
            if len(typmatch) != 1:
                raise TypeError()
            else:
                self.typ = typ
                self.vtyp = typmatch[0][0]
                self.ltyp = typmatch[0][1]
        except TypeError:
            raise io_keyword_error("Invalid typ argument passed to keyword")
        
        if type(description) is not str:
            raise io_keyword_error("Invalid description argument passed to keyword")
        else:
            self.description = description

# Just a convenient wrapper function

def io_freeform_open(fname, keywords=None):
    return io_freeform_file(fname, keywords)

def io_unit_to_atomic(val, units):
    if units not in phys_units:
        raise ValueError("Invalid units " + units)
    return val*phys_units[units]['value']

def io_unit_to_default(val, units):
    if units not in phys_units:
        raise ValueError("Invalid units " + units)
    dim = phys_units[units]['dim']
    return val*phys_units[units]['value']/phys_units[default_units[dim]]['value']

class io_freeform_file(object):
    """IO freeform file class - opens and loads the contents of a file defined by the freeform format.
    
    """
    
    def __init__(self, fname=None, keywords=None):
        """Initialize a io_freeform_file object. Can optionally load a file on initialization

        Args:
        
            fname(str, optional): name of the file to open
            
            keywords([keyword list], optional): a list of instances of the keyword class, defining the keywords that can be found in the file. If not present, all keywords will be accepted
            
        Raises a io_freeform_error if an invalid keyword is present or if blocks are not properly defined.
        """
        
        self.keywords = None
        self.keyvals = collections.OrderedDict()

        if fname is not None:
            self.freeform_load(fname, keywords)

    def __copy__(self):
        new_copy = type(self)()
        new_copy.keywords = copy.copy(self.keywords)
        new_copy.keyvals  = copy.copy(self.keyvals)
        return new_copy 
    
    def freeform_load(self, fname, keywords=None):
        """Open a file and parse its keywords
        
        Args:
        
            fname(str): name of the file to open
            
            keywords([keyword list], optional): a list of instances of the keyword class, defining the keywords that can be found in the file. If not present, all keywords will be accepted
            
        Raises a io_freeform_error if an invalid keyword is present or if blocks are not properly defined.
        """
        
        # First, check that the keywords argument is valid, and open the file
        
        if type(fname) is not str:
            raise io_freeform_error("Invalid fname argument passed to io_freeform_load")
        else:
            try:
                filelines = open(fname, 'r').readlines()
            except IOError:
                raise io_freeform_error("File " + fname + " passed to io_freeform_load not found")
        
        if keywords is not None:
            if type(keywords) is not list or not all([type(k) is keyword for k in keywords]):
                raise io_freeform_error("Invalid keywords argument passed to io_freeform_load")
            else:
                self.keywords = {}
                for k in keywords:
                    self.keywords[k.key.upper()] = k
        else:
            self.keywords = None
        
        # Second, read the file properly and store the contents in an ordered dict to keep the shape of the file
        
        keyw = None
        read_block = False
        
        for i, l in enumerate(filelines):
            
            # Strip all comments, aka anything after a hash
            try:
                i_c = l.index('#')
                l = l[:i_c]
            except ValueError:
                pass    # No comments in this line
            try:
                i_c = l.index('!')
                l = l[:i_c]
            except ValueError:
                pass    # No comments in this line
            
            lsplit = l.replace(':', ' ').split()

            if len(lsplit) == 0:
                # Empty
                continue
            
            
            # Are we reading a block or not?
            
            if read_block:
                if lsplit[0].upper() == '%ENDBLOCK':
                    if len(lsplit) == 1 or lsplit[1].upper() != keyw:
                        raise io_freeform_error("Out of place end of block at line %i in io_freeform_file" % i)
                    else:
                        read_block = False
                else:
                    self.keyvals[keyw] += [' '.join(lsplit)]
            else:
                # Check the first word
                
                # Is it a block?
                read_block = (lsplit[0].upper() == '%BLOCK')
                if read_block:
                    if len(lsplit) == 1:
                        raise io_freeform_error("Unrecognizable block at line %i in io_freeform_file" % i)
                    else:
                        keyw = lsplit[1].upper()
                else:
                    keyw = lsplit[0].upper()
                
                # Is the keyword admissible?
                if self.keywords is not None:
                    if keyw not in self.keywords:
                        raise io_freeform_error("Unrecognizable keyword " + keyw + " at line %i in io_freeform_file" % i)
                    # Is it, correctly, a block?
                    if read_block and self.keywords[keyw].vtyp != 'B':
                        raise io_freeform_error("Unrecognizable block " + keyw + " at line %i in io_freeform_file" % i)

                # Is the keyword duplicated?

                if keyw in self.keyvals:
                    raise io_freeform_error("Duplicated keyword " + keyw + " at line %i in io_freeform_file" % i)
                
                # Now save the value
                
                if read_block:
                    self.keyvals[keyw] = []
                else:
                    self.keyvals[keyw] = ' '.join(lsplit[1:])
            
    
    def freeform_present(self, key):
        """Check if a key is present in the current file
        
        Args:
            
            key(str): the key to check
        
        Returns:
            
            present(bool): True if the key is present
        """
        
        return key.upper() in self.keyvals
    
    def freeform_string(self, key, value=None):
        """Returns or sets the value of a keyword as string. Performs a type check if keywords are present
        
        Args:
            
            key(str): the keyword whose value has to be retrieved
            
            value(str): if present, sets the value of the keyword
        
        Returns:
            
            value(str): the value of the keyword as string
        """
        
        if self.keywords is not None:
            try:
                assert self.keywords[key.upper()].vtyp == 'S'
            except (KeyError, AssertionError):
                raise io_freeform_error("Keyword " + key + " in io_freeform_file is not defined as type String")

        if value is None:
            if not self.freeform_present(key):
                raise io_freeform_error("Keyword " + key + " not present in io_freeform_file")
        else:
            self.keyvals[key.upper()] = str(value)
            
        return self.keyvals[key.upper()]
    
    def freeform_integer(self, key, value=None):
        """Returns or sets the value of a keyword as integer. Performs a type check if keywords are present
        
        Args:
            
            key(str): the keyword whose value has to be retrieved or set
            
            value(int, optional): if present, sets the value of the keyword
        
        Returns:
            
            value(int): the value of the keyword as integer
        """

        if self.keywords is not None:
            try:
                assert self.keywords[key.upper()].vtyp == 'I'
            except (KeyError, AssertionError):
                raise io_freeform_error("Keyword " + key + " in io_freeform_file is not defined as type Integer")
        
        if value is None:
            if not self.freeform_present(key):
                raise io_freeform_error("Keyword " + key + " not present in io_freeform_file")
        else:
            self.keyvals[key.upper()] = str(value)
        
        try:
            return int(self.keyvals[key.upper()])
        except ValueError:
            raise io_freeform_error("Invalid arguments for keyword " + key + " in io_freeform_file")
    
    def freeform_real(self, key, value=None):
        """Returns or sets the value of a keyword as float (real). Performs a type check if keywords are present
        
        Args:
            
            key(str): the keyword whose value has to be retrieved or set
            
            value(float, optional): if present, sets the value of the keyword
        
        Returns:
            
            value(float): the value of the keyword as float
        """
        
        if self.keywords is not None:
            try:
                assert self.keywords[key.upper()].vtyp == 'R'
            except (KeyError, AssertionError):
                raise io_freeform_error("Keyword " + key + " in io_freeform_file is not defined as type Real")

        if value is None:
            if not self.freeform_present(key):
                raise io_freeform_error("Keyword " + key + " not present in io_freeform_file")
        else:
            self.keyvals[key.upper()] = str(value)
        
        try:
            return float(self.keyvals[key.upper()])
        except ValueError:
            raise io_freeform_error("Invalid arguments for keyword " + key + " in io_freeform_file")
        
    def freeform_physical(self, key, dim, value=None, unit=None, usedef=False):
        """Returns or sets the value of a keyword as physical quantity in atomic/default units. Performs a type check if keywords are present
        
        Args:
            
            key(str): the keyword whose value has to be retrieved or set
            
            dim(str): the type of dimension of the keyword
        
            value(int, optional): if present, sets the value of the keyword

            unit(str, optional): if present, sets the units in which value is expressed (otherwise default units are considered)

            usedef(bool, optional): if present and set to True, default units instead of atomic ones are used in output
            
        Returns:
            
            value(float): the value of the keyword as real, in atomic units
        """
        
        if self.keywords is not None:
            try:
                assert self.keywords[key.upper()].vtyp == 'P'
            except (KeyError, AssertionError):
                raise io_freeform_error("Keyword " + key + " in io_freeform_file is not defined as type Physical")
        
        if dim not in default_units:
            raise io_freeform_error("Invalid dimension " + dim + " passed to freeform_physical")

        if unit is not None and unit not in phys_units:
            raise io_freeform_error("Invalid unit " + unit + " passed to freeform_physical")

        if value is None:
            if not self.freeform_present(key):
                raise io_freeform_error("Keyword " + key + " not present in io_freeform_file")
            
            rawvals = self.keyvals[key.upper()].split()
            
            if len(rawvals) > 2:
                raise io_freeform_error("Invalid arguments for keyword " + key + " in io_freeform_file")
            
            try:
                value = float(rawvals[0])
                unit = rawvals[1] if len(rawvals) > 1 else default_units[dim]
            except ValueError:
                raise io_freeform_error("Invalid arguments for keyword " + key + " in io_freeform_file")
        else:
            if unit is None:
                unit = default_units[dim]
            self.keyvals[key.upper()] = str(value) + ' ' + phys_units[unit]['unit']
            
        try:
            if not usedef:
                return io_unit_to_atomic(value, unit)
            else:
                return io_unit_to_default(value, unit)
        except ValueError:
            raise io_freeform_error("Invalid arguments for keyword " + key + " in io_freeform_file")
        
    def freeform_defined(self, key, value=None):
        """Returns or sets the value of whether a keyword is defined or not. Performs a type check if keywords are present
        
        Args:
            
            key(str): the keyword whose value has to be retrieved

            value(bool, optional): if present, sets the value of the keyword
        
        Returns:
            
            value(bool): whether the keyword is defined or not
        """
        
        if self.keywords is not None:
            try:
                assert self.keywords[key.upper()].vtyp == 'D'
            except (KeyError, AssertionError):
                raise io_freeform_error("Keyword " + key + " in io_freeform_file is not defined as type Defined")

        if value is not None:
            if value:
                self.keyvals[key.upper()] = ''
            else:
                try:
                    del(self.keyvals[key.upper()])
                except KeyError:
                    pass    # It never existed to begin with...
        
        return self.freeform_present(key)
            
    def freeform_boolean(self, key, value=None):
        """Returns or sets the value of a keyword as boolean (logical). Performs a type check if keywords are present
        
        Args:
            
            key(str): the keyword whose value has to be retrieved
        
            value(bool, optional): if present, sets the value of the keyword
            
        Returns:
            
            value(bool): the value of the keyword as bool
        """
        
        if self.keywords is not None:
            try:
                assert self.keywords[key.upper()].vtyp == 'L'
            except (KeyError, AssertionError):
                raise io_freeform_error("Keyword " + key + " in io_freeform_file is not defined as type Logical")

        if value is None:
            if not self.freeform_present(key):
                raise io_freeform_error("Keyword " + key + " not present in io_freeform_file")
        else:
            self.keyvals[key.upper()] = {True: 'TRUE', False: 'FALSE'}[value]
            
        try:
            return {'TRUE': True, 'FALSE': False}[self.keyvals[key.upper()].strip().upper()]
        except KeyError:
            raise io_freeform_error("Invalid arguments for keyword " + key + " in io_freeform_file")
            
    def freeform_real_vector(self, key, value=None):
        """Returns the value of a keyword as list of floats (real vector). Performs a type check if keywords are present
        
        Args:
            
            key(str): the keyword whose value has to be retrieved
        
            value(list, optional): if present, sets the value of the keyword
            
        Returns:
            
            value(list): the value of the keyword as a list
        """
        
        if self.keywords is not None:
            try:
                assert self.keywords[key.upper()].vtyp == 'V'
            except (KeyError, AssertionError):
                raise io_freeform_error("Keyword " + key + " in io_freeform_file is not defined as type Real Vector")

        if value is None:
            if not self.freeform_present(key):
                raise io_freeform_error("Keyword " + key + " not present in io_freeform_file")
        else:
            try:
                self.keyvals[key.upper()] = '%f %f %f' % value
            except TypeError:
                raise ValueError("value must be a list/tuple of 3 elements")
            
        try:
            val = [float(x) for x in self.keyvals[key.upper()].split()]
            if len(val) != 3:
                raise ValueError()
            return val
        except ValueError:
            raise io_freeform_error("Invalid arguments for keyword " + key + " in io_freeform_file")
        
    def freeform_integer_vector(self, key, value=None):
        """Returns or sets the value of a keyword as list of integers (integer vector). Performs a type check if keywords are present
        
        Args:
            
            key(str): the keyword whose value has to be retrieved
        
            value(list, optional): if present, sets the value of the keyword
            
        Returns:
            
            value(list): the value of the keyword as a list
        """
        
        if self.keywords is not None:
            try:
                assert self.keywords[key.upper()].vtyp == 'W'
            except (KeyError, AssertionError):
                raise io_freeform_error("Keyword " + key + " in io_freeform_file is not defined as type Integer Vector")

        if value is None:
            if not self.freeform_present(key):
                raise io_freeform_error("Keyword " + key + " not present in io_freeform_file")
        else:
            try:
                self.keyvals[key.upper()] = '%i %i %i' % value
            except TypeError:
                raise ValueError("value must be a list/tuple of 3 elements")
                
            
        try:
            val = [int(x) for x in self.keyvals[key.upper()].split()]
            if len(val) != 3:
                raise ValueError()
            return val
        except ValueError:
            raise io_freeform_error("Invalid arguments for keyword " + key + " in io_freeform_file")

    def freeform_block(self, key, value=None):
        """Returns or sets the value of a keyword as list of strings (block data, one per line). Performs a type check if keywords are present
        
        Args:
            
            key(str): the keyword whose value has to be retrieved

            value(list, optional): if present, sets the value of the keyword
        
        Returns:
            
            value(list): the value of the data block as a list of strings
        """
        
        if self.keywords is not None:
            try:
                assert self.keywords[key.upper()].vtyp == 'B'
            except (KeyError, AssertionError):
                raise io_freeform_error("Keyword " + key + " in io_freeform_file is not defined as type Block")

        if value is None:
            if not self.freeform_present(key):
                raise io_freeform_error("Keyword " + key + " not present in io_freeform_file")
        else:
            try:
                self.keyvals[key.upper()] = [str(l) for l in value]
            except TypeError:
                raise ValueError("value must be a list/tuple")
            
        return self.keyvals[key.upper()]

    def freeform_remove(self, key):
        """Removes an existing keyword

        Args:

            key(str): the keyword that needs to be removed

        Returns:

            None
        """

        if self.keywords is not None:
            if key.upper() not in self.keywords:
                raise io_freeform_error("Keyword " + key + " is not defined for io_freeform_file")

        try:
            del(self.keyvals[key.upper()])
        except KeyError:
            pass

    def freeform_print(self):
        """Prints out a freeform file

        Args:

            None

        Returns:

            printed_file(str): the printed out freeform file, as string
        """

        printed_file = ""

        for k in self.keyvals:

            # Check if it's a simple keyword or a block

            kval = self.keyvals[k]

            if type(kval) == list:    # It's a block

                printed_file += "%BLOCK {0}\n".format(k)
                for l in kval:
                    printed_file += "%s\n" % l
                printed_file += "%ENDBLOCK {0}\n".format(k)
            else:
                printed_file += "%s\t%s\n" % (k, kval)

            printed_file += "\n"

        return printed_file

if __name__ == "__main__":

    f = io_freeform_file()

    print f.freeform_print()