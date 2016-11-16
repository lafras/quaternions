from collections import namedtuple

import decimal
Decimal = decimal.Decimal

_PREC_DEC = 2**5
_PREC_DEC_EQ = 2**4
_PREC_DEC_REPR = 2

decimal.setcontext(decimal.Context(prec=_PREC_DEC))
context = decimal.getcontext()

import fractions
math = fractions.math
operator = fractions.operator


PreQuaternion = namedtuple('H','a b c d',verbose=False)

class Quaternion(PreQuaternion):
    
    def __new__(_cls, a, b, c, d):
        a,b,c,d = map(context.create_decimal, map(str,[a,b,c,d]))
        return _cls.__bases__[0].__new__(_cls, a,b,c,d)

    def __repr__(self):
        return "{a:-.2f}{b:+.2f}i{c:+.2f}j{d:+.2f}k".format(**self._asdict())

    def __abs__(self):
        return self.norm()

    def __add__(self,other):
        a1,b1,c1,d1 = self
        a2,b2,c2,d2 = other
        return Quaternion(a1+a2,b1+b2,c1+c2,d1+d2)

    def __cmp__(self,other):
        raise TypeError("No ordering relation defined on quaternions")
      
    def __truediv__(self,other):
                
        if isinstance(other, Quaternion):
            return operator.mul(self, other.inv())
        else:
            try:
                f = Decimal(str(other))
            except:                
                raise NotImplementedError
            a,b,c,d = self
            return Quaternion(a/f,b/f,c/f,d/f)

    def __rdiv__(self,other):
        raise NotImplementedError

    def __eq__(self,other):
        
        with decimal.localcontext() as ctx:
            ctx.prec = _PREC_DEC_EQ
            eq = all(ctx.compare(a,b).is_zero() for a,b in zip(self,other))
            
        return eq

    def __mul__(self,other):
        
        if isinstance(other, Quaternion):
            a1,b1,c1,d1 = self
            a2,b2,c2,d2 = other
            return Quaternion(a1*a2 - b1*b2 - c1*c2 - d1*d2,
                              a1*b2 + b1*a2 + c1*d2 - d1*c2,
                              a1*c2 - b1*d2 + c1*a2 + d1*b2,
                              a1*d2 + b1*c2 - c1*b2 + d1*a2)
        else:
            try:
                f = Decimal(str(other))
            except:
                raise NotImplementedError
            a,b,c,d = self
            return Quaternion(f*a,f*b,f*c,f*d)    
    
    def __rmul__(self,other):
        return operator.mul(self,other)

    def __neg__(self):
        a,b,c,d = self
        return Quaternion(-a,-b,-c,-d)
    
    def __sub__(self,other):
        a1,b1,c1,d1 = self
        a2,b2,c2,d2 = other
        return Quaternion(a1-a2,b1-b2,c1-c2,d1-d2)

    def conj(self):
        a,b,c,d = self
        return Quaternion(a,-b,-c,-d)
        
    def inv(self):
        a,b,c,d = self.conj()
        n = self.norm() ** 2
        return Quaternion(a/n,b/n,c/n,d/n)
    
    def matrix(self):
        """ [[ a, b, c, d],
             [-b, a,-d, c],
             [-c, d, a,-b],
             [-d,-c, b, a]] """
             
        raise NotImplementedError
    
    def norm(self):
        a,b,c,d = self
        return sum([a**2,b**2,c**2,d**2]).sqrt()
    
    def normalise(self):
        return self / self.norm()

    def polar(self):
        raise NotImplementedError

    def scalar(self):
    
        """The scalar part as a Quaternion"""
    
        return Quaternion(self.a,0.0,0.0,0.0)
    
    def vector(self):
        
        """The vector part as a Quaternion"""
        
        return Quaternion(0.0,self.b,self.c,self.d)
    
    def unit(self):
        
        """Return a unit vector pointing in the direction of vec(q)"""
        
        v = self.vector()
        return v.norm()

def cross(p,q):
    
    """The cross product of two purely imaginary quaternions"""
    
    try:
        assert (p.a == q.a == 0)
    except AssertionError:
        raise ValueError("Arguments are not purely imaginary")
    
    _,b1,c1,d1 = p
    _,b2,c2,d2 = q
    
    return Quaternion(0,c1*d2-d1*c2,d1*b2-b1*d2,b1*c2-c1*b2)

def dot(p,q):
    
    """The dot product of two purely imaginary quaternions"""
    
    try:
        assert (p.a == q.a == 0)
    except AssertionError:
        raise ValueError("Arguments are not purely imaginary")
    
    _,b1,c1,d1 = p
    _,b2,c2,d2 = q
    
    return Quaternion(b1*b2+c1*c2+d1*d2,0,0,0)
    
def exp(q):

    """The exponential of a quaternion"""
    
    v = q.vector()
    w = v / v.norm() * math.sin(v.norm())
    s = Quaternion(math.cos(v.norm()),0.0,0.0,0.0)        

    return math.exp(q.a) * (s+w)

def log(q):
    
    """The natural logarithm of a quaternion"""    

    v = q.vector()
    w = v / v.norm() * math.acos(q.a/q.norm())
    s = Quaternion(math.log(q.norm()),0.0,0.0,0.0)
    
    return s+w
        
    
if __name__ == '__main__':
    
#    p = Quaternion(1.0,0.1,0.2,0.3)
    
    import random
    
    p = Quaternion(*(random.uniform(-1,1) for i in range(4)))
    q = Quaternion(*(random.uniform(-1,1) for i in range(4))) 
    
    print(p)
#    print(q)
    

    print(log(exp(p)))
    print(exp(log(p)))
    
    print(log(exp(p)) == exp(log(p)))
