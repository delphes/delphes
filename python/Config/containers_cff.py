##
# \file
# \ingroup python
# \note Shamelessly self-stolen from CepGen/TDAnalyser
# \author Laurent Forthomme <laurent.forthomme@cern.ch>

class PrintHelper(object):
    """Helper class for the pretty-printing of configuration parameters"""
    _indent = 0
    _indent_size = 4

    def indent(self):
        """Move to the next indentation block"""
        self._indent += self._indent_size

    def unindent(self):
        """Go up to the previous indentation block"""
        self._indent -= self._indent_size

    def indentation(self):
        """Current indentation level"""
        return ' '*self._indent


class Parameters(dict):
    """A raw list of steering parameters"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __init__(self, *args, **kwargs):
        """Construct a parameters set from dictionary arguments

        Args:
            args: list of arguments

        Kwargs:
            kwargs: list of named (keyworded) arguments

        Examples:
            >>> print(dict(Parameters(foo = 'bar', foo2 = 42)))
            {'foo': 'bar', 'foo2': 42}
        """
        self.update(*args, **kwargs)
        super(Parameters, self).__init__(*args, **kwargs)

    def __deepcopy__(self, memo):
        """Override the default dict deep copy operator"""
        from copy import deepcopy
        return Parameters([(deepcopy(k, memo), deepcopy(v, memo)) for k, v in self.items()])

    def dump(self, printer=PrintHelper()):
        """Human-readable dump of this object"""
        out = self.__class__.__name__+'(\n'
        printer.indent()
        for k, v in self.items():
            out += ('%s%s = ' % (printer.indentation(), k))
            if v.__class__.__name__ in ['Parameters', 'Module']:
                out += v.dump(printer)
            elif v.__class__.__name__ in ['list', 'tuple']:
                out += v.__class__.__name__ + '('
                printer.indent()
                for it in v:
                    if it.__class__.__name__ in ['Parameters', 'Module']:
                        out += '\n' + printer.indentation()
                        out += it.dump(printer)
                    else:
                        out += it.__repr__()
                    if it != v[-1]:
                        out += ', '
                out += ')'
                printer.unindent()
            else:
                out += v.__repr__()
            out += ',\n'
        printer.unindent()
        out += printer.indentation()+')'
        return out

    def __repr__(self):
        """Human-readable version of this object"""
        return self.dump()

    def clone(self, *args, **kwargs):
        """Return a deep copy of this object"""
        from copy import deepcopy
        out = deepcopy(self)
        for k in kwargs:
            out[k] = kwargs.get(k)
        return type(self)(out)

    def load(self, mod):
        """Extend this object by an include"""
        from sys import modules
        mod = mod.replace('/', '.')
        __import__(mod)
        self.extend(modules[mod])


class Module(Parameters):
    """A named parameters set to steer a generic module

    Attributes:
        type: Name of this module
    """

    def __init__(self, type_, *args, **kwargs):
        """Construct a module parameters set from dictionary arguments

        Args:
            type_: module type
            args: list of arguments

        Kwargs:
            kwargs: list of named (keyworded) arguments

        Examples:
            >>> print(dict(Module('module1', foo = 'bar', foo2 = 42)))
            {'foo': 'bar', 'foo2': 42, 'type': 'module1'}
        """
        super(Module, self).__init__(*args, **kwargs)
        self.type = type_

    def __len__(self):
        """Number of keys handled"""
        return dict.__len__(self) - 1  # discard the type key

    def dump(self, printer=PrintHelper()):
        """Human-readable dump of this object"""
        out = self.__class__.__name__+'('+self.type.__repr__()+',\n'
        mod_repr = self.clone('')
        mod_repr.pop('type', None)
        for k, v in mod_repr.items():
            printer.indent()
            out += ('%s%s = ' % (printer.indentation(), k))
            if v.__class__.__name__ not in ['Parameters', 'Module']:
                out += v.__repr__()
            else:
                out += v.dump(printer)
            out += ',\n'
            printer.unindent()
        out += printer.indentation()+')'
        return out

    def __repr__(self):
        """Human-readable version of this object"""
        return self.dump()

    def clone(self, type_='', **kwargs):
        """Return a deep copy of this object"""
        out = Parameters(self).clone(**kwargs)
        if type_:
            out.type = type_
        return out


class Sequence(list):
    """An ordered modules sequence"""
    MODULE = object()

    def __init__(self, *args):
        """Construct an ordered sequence of modules from a list

        Args:
            args: list of modules

        Examples:
            >>> module1 = Module('test1', foo = 'bar')
            >>> module2 = Module('test2', foo2 = 42)
            >>> print(Sequence(module1, module2))
            Sequence([Module('test1',
                foo = 'bar',
            ), Module('test2',
                foo2 = 42,
            )])
        """
        super(Sequence, self).__init__(args)

    def __delitem__(self, index):
        """Remove an element from the sequence"""
        self[index] = self.MODULE

    def __iter__(self):
        """Iterator definition for the sequence"""
        return (item for item in super().__iter__() if item is not self.MODULE)

    def __eq__(self, other):
        """Equality operator"""
        if isinstance(other, Sequence):
            return all(x == y for x, y in zip(self, other))
        return super().__eq__(other)

    def __repr__(self):
        """Human-readable representation of this sequence"""
        return type(self).__name__+'('+super(Sequence, self).__repr__()+')'


if __name__ == '__main__':
    import unittest

    class TestTypes(unittest.TestCase):
        """Small collection of tests for our new types"""

        def testModules(self):
            """Test the Module object"""
            mod = Module('empty')
            self.assertEqual(len(mod), 0)
            mod.param1 = 'foo'
            self.assertEqual(len(mod), 1)
            # playing with module clones
            mod_copy = mod.clone('notEmpty', param1 = 'boo', param2 = 'bar')
            self.assertEqual(mod.param1, 'foo')
            self.assertEqual(mod_copy.param1, 'boo')
            self.assertEqual(mod_copy.param2, 'bar')
            self.assertEqual(mod.param1+mod_copy.param2, 'foobar')

        def testParameters(self):
            """Test the Parameters object"""
            params = Parameters(
                first = 'foo',
                second = 'bar',
                third = 42,
                fourth = (1, 2),
            )
            params_copy = params.clone(
                second = 'bak',
            )
            self.assertEqual(len(params), 4)
            self.assertEqual(params.first, params['first'])
            self.assertEqual(params['second'], 'bar')
            self.assertTrue(int(params.third) == params.third)
            self.assertEqual(len(params.fourth), 2)
            self.assertEqual(params.second, 'bar')
            # playing with parameter clones
            self.assertEqual(params_copy.second, 'bak')
            # check that the clone does not change value if the origin does
            # (i.e. we indeed have a deep copy and not a shallow one...)
            params.third = 43
            self.assertEqual(params.third, 43)
            self.assertEqual(params_copy.third, 42)

        def testSequences(self):
            """Test the Sequence object"""
            mod1 = Module('module1',
                foo = 'bar',
                foo2 = 42
            )
            mod2 = Module('module2',
                foo = [i for i in range(4)]
            )
            seq = Sequence(mod1, mod2)
            self.assertEqual(len(seq), 2)
            # check that all sequence modules have their attributes
            # conserved in the sequence construction
            self.assertEqual(len(seq[0]), 2)
            self.assertEqual(seq[0].foo, 'bar')
            self.assertEqual(seq[0].type, 'module1')
            self.assertEqual(len(seq[1]), 1)
            self.assertEqual(seq[1].type, 'module2')
            self.assertEqual(len(seq[1].foo), 4)

    unittest.main()
