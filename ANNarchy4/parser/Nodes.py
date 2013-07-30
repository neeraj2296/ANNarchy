from Definitions import *
import Tree
            
## Nodes
class TreeBuilder:

    def __init__(self, machine, expr_list):
        self.machine = machine
        self.expr_list = expr_list
    
    def build(self):
        # Split the list around the equal sign
        id_equal = self.expr_list.index(EQUAL)
        self.left_group, self.right_group = Group(self.machine, self.expr_list[:id_equal]), Group(self.machine, self.expr_list[id_equal+1:])
        # Start the group analysis
        for group in [self.left_group, self.right_group]:
            group.group()
            group.analyse()       
        # Return ordered tree
        tree = Equal(self.machine, left=self.left_group, right=self.right_group).build_hierarchy()
        return tree
        
# Node class
class Node:
    def __init__(self, machine):
        self.machine = machine
        
    def __repr__(self):
        return self.cpp()

    def hierarchize(self, side):
        """ Create a hierarchy using the recognized nodes. """
        if len(side)==1: # The group has only one item
            if isinstance(side[0], Group): # A term between brackets.
                newside = self.hierarchize(side[0])
                newside.parent= self
                return newside
            else: # Terminal leaf
                newside = side[0]
                newside.parent= self
                return newside
        found = False
        for i in range(len(side)): # if then else
            if isinstance(side[i], If):
                newside = If(self.machine, value = 'if', child=(side[i+1:]) )
                newside.parent= self
                return newside
        for i in range(len(side)): # logical operator
            if isinstance(side[i], Logical):
                newside = Logical(self.machine, 
                                    value=side[i].value,
                                    left=(side[:i]), 
                                    right=(side[i+1:]) )
                newside.parent= self
                return newside
        for i in range(len(side)): # not operator
            if isinstance(side[i], Not):
                newside = Not(self.machine,
                                value = side[i].value,
                                child = side[i+1:] )
                newside.parent= self
                return newside
        for i in range(len(side)): # comparator
            if isinstance(side[i], Comparator):
                newside = Comparator(self.machine, 
                                    value=side[i].value,
                                    left=(side[:i]), 
                                    right=(side[i+1:]) )
                newside.parent= self
                return newside
        for i in range(len(side)): # add
            if isinstance(side[i], PlusOperator):
                newside = PlusOperator(self.machine, 
                                    left=(side[:i]), 
                                    right=(side[i+1:]) )
                newside.parent= self
                found=True
                break
        if not found:
            for i in range(len(side)): # sub
                if isinstance(side[i], MinusOperator):
                    if i >0:
                        newside = MinusOperator(self.machine, 
                                            left=(side[:i]), 
                                            right=(side[i+1:]) )
                    else:
                        newside = MinusOperator(self.machine, 
                                            left=[Constant(self.machine, value='')], 
                                            right=(side[i+1:]) )
                    newside.parent= self
                    found=True
                    break
        if not found: # Multiply
            for i in range(len(side)):
                if isinstance(side[i], MultOperator):
                    newside = MultOperator(self.machine, 
                                        left=(side[:i]), 
                                        right=(side[i+1:]) )
                    newside.parent= self
                    found=True
                    break
        if not found: # Divide
            for i in range(len(side)-1, 0, -1):
                if isinstance(side[i], DivOperator):
                    newside = DivOperator(self.machine, 
                                        left=(side[:i]), 
                                        right=(side[i+1:]) )
                    newside.parent= self
                    found=True
                    break
        if not found: # power function
            for i in range(len(side)):
                if isinstance(side[i], PowOperator):
                    newside = PowOperator(self.machine, 
                                        left=(side[:i]), 
                                        right=(side[i+1:]) )
                    newside.parent= self
                    found=True
                    break
        if not found: # No operators, only functions are left
            for i in range(len(side)):
                if isinstance(side[i], Function):
                    newside = Function(self.machine,
                                    value = side[i].value,
                                    child = side[i+1] )
                    newside.parent= self
                    found=True
                    break
        if not found:
            print self.machine.expr
            print 'Error: could not analyse the expression'
            exit(0)
        return newside
        
# Equal object        
class Equal(Node):
    def __init__(self, machine, left, right):
        Node.__init__(self, machine)
        self.machine = machine
        self.left = left
        self.right = right
        self.order=0
                
    def cpp(self):
        if self.order == 1:
            return self.left.cpp() + " += dt_*" + self.right.cpp()
        return self.left.cpp() + " = " + self.right.cpp()
                
    def latex(self):
        return self.left.latex() + " = " + self.right.latex()
        
    def build_hierarchy(self):
        left=None
        right=None
        # left member
        if len(self.left) == 1:
            if not isinstance(self.left[0], MainVariable) and not isinstance(self.left[0], Gradient):
                print self.machine.expr 
                print 'Error: the left term should be the updated variable'
                exit(0)
            else:
                left = self.left[0]
        else: # there is at least one operator
            left = self.hierarchize(self.left) 
        self.left = left
        self.left.parent= self  
                    
        # right member 
        right = self.hierarchize(self.right) 
        self.right = right
        self.right.parent= self
        return self
        
# Leaf (no need to go further)
class Leaf:
    def __repr__(self):
        return self.cpp()       
        
# Constants
class Constant(Leaf):
    def __init__(self, machine, value):
        self.machine = machine
        self.value = value
        
    def cpp(self):
        if self.value in PI:
            return ' ' + cpp_equivalent(self.value) + ' '
        elif self.value in TRUE:
            return ' ' + cpp_equivalent(self.value) + ' '
        elif self.value in FALSE:
            return ' ' + cpp_equivalent(self.value) + ' '
        elif self.value == '': # case of negative numbers
            return ''
        return ' ' + str(float(self.value)) + ' '
        
    def latex(self):
        if self.value in PI:
            return '{\\pi}'
        return '{' + str(self.value) + '}'
        
class Variable(Leaf):
    def __init__(self, machine, value):
        self.machine = machine
        self.value = value
    def cpp(self):
        if self.value == 't':
            return ' time_ '
        return ' ' + str(self.value)+'[i] '
    def latex(self):
        return '{\\text{'+str(self.value) + '}}'
        
class MainVariable(Variable):
    def __init__(self, machine, value):
        Variable.__init__(self, machine, value)
        
class Parameter(Leaf):
    def __init__(self, machine, value):
        self.machine = machine
        self.value = value
    def cpp(self):
        return str(self.value)
    def latex(self):
        return '{\\text{'+str(self.value) + '}}'
        
class Gradient(Leaf):
    def __init__(self, machine, value):
        self.machine = machine
        self.value = value
    def cpp(self):
        return ' ' + str(self.value)+'[i] '
    def latex(self):
        return '{\\frac{d\\text{'+str(self.value) + '}}{dt}}'
        
class PSP(Leaf):
    def __init__(self, machine, value):
        self.machine = machine
        self.value = value[0]
    def cpp(self):
        if not self.value in self.machine.targets:
            print self.machine.expr
            print 'Error: the target', self.value, 'does not exist on this neuron. The sum will be 0.0.'
            return ' 0.0 '
        return ' sum(i, ' + str(self.machine.targets.index(self.value))+') '
    def latex(self):
        return '{\\sum_{i}^{'+str(self.value)+'} w_{i} \\cdot \\text{rate}_{i}}'
        
class PreVariable(Leaf):
    def __init__(self, machine, value):
        self.machine = machine
        self.value = value
    def cpp(self):
        variable =  self.value.split('.')[1]
        if variable=='rate':
            return ' (*pre_rates_) [ rank_ [i] ] '
        return ' pre_population_->get'+variable.capitalize()+'( rank_[i] ) '
    def latex(self):
        return '{\\text{'+str(self.value)+'}}'
        
class PostVariable(Leaf):
    def __init__(self, machine, value):
        self.machine = machine
        self.value = value
    def cpp(self):
        variable =  self.value.split('.')[1]
        if variable=='rate':
            return ' (*post_rates_) [ post_neuron_rank_ ] '
        return ' post_population_->get'+variable.capitalize()+'( post_neuron_rank_ ) '
    def latex(self):
        return '{\\text{'+str(self.value)+'}}'


        
# Unitary operator    
class Function(Node):
    def __init__(self, machine, value, child=None, hierarchize=True):
        Node.__init__(self, machine)
        self.machine = machine
        self.value = value
        self.child = child
        if self.child != None and hierarchize:
            self.child = self.hierarchize(self.child)
            self.child.parent=self
            
    def cpp(self):
        if self.child != None:
            return str(cpp_equivalent(self.value))+'('+self.child.cpp()+')'  
        return str(self.value)+'()'
        
    def latex(self):
        if self.child != None:
            if self.value in POS:
                return '{('+self.child.latex()+')^+}' 
            elif self.value in NEG:
                return '{('+self.child.latex()+')^-}' 
            elif self.value in SQRT:
                return '{' + latex_equivalent(self.value)+'{'+self.child.latex()+'}}' 
            else:
                return '{' + latex_equivalent(self.value)+'{('+self.child.latex()+')}}'  
        return '{\\text{'+str(self.value)+'}()}'
        
# Not operator    
class Not(Node):
    def __init__(self, machine, value, child=None, hierarchize=True):
        Node.__init__(self, machine)
        self.machine = machine
        self.value = value
        self.child = child
        print self.child
        if self.child != None and hierarchize:
            self.child = self.hierarchize(self.child)
            self.child.parent=self
            
    def cpp(self):
        if self.child != None:
            return str(cpp_equivalent(self.value))+'('+self.child.cpp()+')'  
        return str(self.value)+'()'
        
    def latex(self):
        return '{\\bar{'+str(self.value)+'}}'
        
# Operator
class Operator(Node):
    def __init__(self, machine, value, left, right, hierarchize=True):
        Node.__init__(self, machine)
        self.machine = machine
        self.value = value
        self.left= left
        self.right=right
        if self.left != None and  hierarchize:
            self.left = self.hierarchize(self.left)
            self.left.parent = self
        if self.right != None and hierarchize:
            self.right = self.hierarchize(self.right)
            self.right.parent = self
    
    def cpp(self):
        if self.left != None and self.right != None:
            return ' (' + self.left.cpp() + str(self.value) + self.right.cpp() + ') '
        return ' (' + str(self.value) + ') '

# Ternary operator        
class If(Node):
    def __init__(self, machine, value, child=None, hierarchize=True):
        self.machine = machine
        self.value = value
        self.child = child
        if self.child is not None and hierarchize:
            id_then = -1
            id_else = -1
            nb_if = 0
            for rk, item in enumerate(self.child):
                if isinstance(item, If):
                    nb_if += 1
                elif isinstance(item, Then):
                    if nb_if == 0:
                        id_then =  rk
                elif isinstance(item, Else):                    
                    if nb_if == 0:
                        id_else =  rk
                    else:
                        nb_if -= 1
                    
            if id_then is -1 or id_else is -1:
                print 'Error in analysing', self.machine.expr
                print 'The conditional should use the (if A then B else C) structure.'
                exit(0)
            self.cond = self.hierarchize(Group(self.machine, self.child[:id_then]))
            self.then = self.hierarchize(Group(self.machine, self.child[id_then+1: id_else]))
            self.elsestatement = self.hierarchize(Group(self.machine, self.child[id_else+1:]))
            
    def cpp(self):
        if self.child:
            return '( (' +  self.cond.cpp() + ') ? ' + self.then.cpp() + ' : ' + self.elsestatement.cpp() + ')'
        else: 
            return ' if ' 
    def latex(self):
        return ''      
        
class Then(Node):
    def __init__(self, machine, value, child=None, hierarchize=True):
        self.machine = machine
        self.value = value
        self.child = child
            
    def cpp(self):
        return ' then ' 
    def latex(self):
        return ''   
           
class Else(Node):
    def __init__(self, machine, value, child=None, hierarchize=True):
        self.machine = machine
        self.value = value
        self.child = child
            
    def cpp(self):
        return ' else ' 
    def latex(self):
        return ''          

# Comparators
class Comparator(Operator):
    def __init__(self, machine, value='<', left=None, right=None, hierarchize=True):
        Operator.__init__(self, machine, value, left, right, hierarchize)
        
    def cpp(self): # comparators are already in C++
        if self.left is not None:
            return self.left.cpp() + cpp_equivalent(self.value) + self.right.cpp() 
        return str(self.value)
            
    def latex(self):
        if self.value is '>':
            return '{>}'
        elif self.value is '>=':
            return '{\geq}'
        elif self.value is '<':
            return '{<}'
        elif self.value is '<=':
            return '{\leq}'
        elif self.value is '!=':
            return '{\neq}'
        return '{' + str(self.value) + '}'
        
# Logical operator
class Logical(Operator):
    def __init__(self, machine, value='and', left=None, right=None, hierarchize=True):
        Operator.__init__(self, machine, value, left, right, hierarchize)
        
    def cpp(self): # comparators are already in C++
        if self.left != None:
            return ' (' + self.left.cpp() + ') ' + str(cpp_equivalent(self.value)) + ' (' + self.right.cpp() + ') '
        return str(self.value)
            
    def latex(self):
        return "\\text{and}"

# Binary operators    
class PlusOperator(Operator):   
    def __init__(self, machine, value='+', left=None, right=None, hierarchize=True):
        Operator.__init__(self, machine, value, left, right, hierarchize)
        
    def latex(self):
        if self.left != None and self.right != None:
            if not isinstance(self.parent, MultOperator):
                return ' {' + self.left.latex() + str(self.value) + self.right.latex() + '} '
            else:
                return ' {(' + self.left.latex() + str(self.value) + self.right.latex() + ')} '
        return ''
    
class MinusOperator(Operator):   
    def __init__(self, machine, value='-', left=None, right=None, hierarchize=True):
        Operator.__init__(self, machine, value, left, right, hierarchize)
    def latex(self):
        if self.left != None and self.right != None:
            if not isinstance(self.parent, MultOperator):
                return ' {' + self.left.latex() + str(self.value) + self.right.latex() + '} '
            else:
                return ' {(' + self.left.latex() + str(self.value) + self.right.latex() + ')} '
        return ''
    
class MultOperator(Operator):   
    def __init__(self, machine, value='*', left=None, right=None, hierarchize=True):
        Operator.__init__(self, machine, value, left, right, hierarchize)
    def latex(self):
        if self.left != None and self.right != None:
            return ' {' + self.left.latex() + ' \cdot ' + self.right.latex() + '} '
        return ''
        
class DivOperator(Operator):   
    def __init__(self, machine, value='/', left=None, right=None, hierarchize=True):
        Operator.__init__(self, machine, value, left, right, hierarchize)
    def latex(self):
        if self.left != None and self.right != None:
            return ' {\\frac{' + self.left.latex() + '}{' + self.right.latex() + '}} '
        return ''
        
class PowOperator(Operator):   
    def __init__(self, machine, value='^', left=None, right=None, hierarchize=True):
        Operator.__init__(self, machine, value, left, right, hierarchize)

    def cpp(self):
        return ' pow('+str(self.left) + ' , ' + str(self.right)+') '

    def latex(self):
        return ' {('+self.left.latex() + ')}^{' + self.right.latex()+'} '

        
# Group of expressions with parenthesis                
class Group(Node):

    def __init__(self, machine, items):
        Node.__init__(self, machine)
        self.machine = machine
        self.items = items
        
    def __repr__(self):
        return str(self.items)
        
    def cpp(self):
        return str(self.items)
        
    def latex(self):
        return str(self.items)
        
    def __len__(self):
        return len(self.items)
        
    def __getitem__(self, k):
        return self.items[k]
        
    def group(self):
        complete = False
        while not complete:
            rk = 0
            idx_left = idx_right = 0
            complete=True
            # Find first-level brackets
            for idx in range(len(self.items)):
                it = self.items[idx]
                if it == BRACKET_LEFT:
                    complete=False
                    rk += 1
                    if rk == 1: # first bracket
                        idx_left = idx
                elif it == BRACKET_RIGHT:
                    rk -= 1
                    if rk== 0: # found the last matching bracket
                        idx_right = idx
                        break
            # Expand
            if idx_right != idx_left:
                subgroup = Group(self.machine, self.items[idx_left+1:idx_right])
                items=[]
                for i in range(idx_left):
                    items.append(self.items[i])
                items.append(subgroup)
                for i in range(idx_right+1, len(self.items)):
                    items.append(self.items[i])
                subgroup.group()
                self.items = items
            
    def analyse(self):
        items=[]
        # Turn the items into objects
        iter_items = iter(self.items)
        for item in iter_items:
            found = False
            if isinstance(item, str):
                # Check if it is a constant
                if belongs_to(item, CONSTANTS):
                    items.append(Constant(self.machine, item))
                    found = True
                # Check if it is an integer or float
                if Tree.isFloat(item) or Tree.isInt(item):
                    items.append(Constant(self.machine, item))
                    found=True
                # Check if it is a parameter
                if item in self.machine.parameters:
                    items.append(Parameter(self.machine, item))
                    found=True
                # Check if it is a variable
                if item in self.machine.other_variables:
                    items.append(Variable(self.machine, item))
                    found=True
                # Check if it is a function
                if belongs_to(item, FUNCTIONS):
                    items.append(Function(self.machine, item))
                    found = True
                # Check if it is a ternary operator
                if belongs_to(item, IF):
                    items.append(If(self.machine, item))
                    found = True
                if belongs_to(item, THEN):
                    items.append(Then(self.machine, item))
                    found = True
                if belongs_to(item, ELSE):
                    items.append(Else(self.machine, item))
                    found = True
                # Check if it is a weighted sum
                if belongs_to(item, SUMS):
                    items.append(PSP(self.machine, iter_items.next()))
                    found = True
                # Check if it is a pre variable
                if not item.find('pre.') == -1:
                    items.append(PreVariable(self.machine, item))
                    found = True
                # Check if it is a post variable
                if not item.find('post.') == -1:
                    items.append(PostVariable(self.machine, item))
                    found = True
                # Check if it is the updated variable
                if item == self.machine.name:
                    items.append(MainVariable(self.machine, item))
                    found=True
                if item == 'd' + self.machine.name: # found the first item of a gradient
                    items.append(Gradient(self.machine, self.machine.name))
                    if not iter_items.next() == '/' or not iter_items.next() == 'dt':
                        print 'Error: the gradient on', self.machine.name, 'is badly expressed, it should be', 'd'+self.machine.name+'/dt'
                        exit(0)
                    found=True
                # Check if it is an operator
                if belongs_to(item, PLUS):
                    items.append(PlusOperator(self.machine, item))
                    found=True
                if belongs_to(item, MINUS):
                    items.append(MinusOperator(self.machine, item))
                    found=True
                if belongs_to(item, MULT):
                    items.append(MultOperator(self.machine, item))
                    found=True
                if belongs_to(item, DIV):
                    items.append(DivOperator(self.machine, item))
                    found=True  
                if belongs_to(item, POW):
                    items.append(PowOperator(self.machine, item))
                    found=True  
                # Check if it is a comparator
                if belongs_to(item, COMPARATORS):
                    items.append(Comparator(self.machine, item))
                    found = True  
                # Check if it is a logical operator
                if belongs_to(item, LOGICALS):
                    items.append(Logical(self.machine, item))
                    found = True  
                # Check if it is a not operator
                if belongs_to(item, NOT):
                    items.append(Not(self.machine, item))
                    found = True                
            # Check if it is a subgroup
            elif isinstance(item, Group):
                item.analyse()
                items.append(item)
                found = True          
            # Else: don't know what to do with it...
            if not found:
                print self.machine.expr
                print 'Error: the term', item, 'is neither part of the allowed vocabulary or a parameter of your model.'
                exit(0)
        self.items = items

