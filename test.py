import unittest

class MultiplicationRuleTestCase(unittest.TestCase):
    
    def setUp(self):
    
        import quaternion
        
        self.i = quaternion.Quaternion(0,1,0,0)
        self.j = quaternion.Quaternion(0,0,1,0)
        self.k = quaternion.Quaternion(0,0,0,1)
        
    def test_ijk(self):
        lhs = self.i*self.j
        rhs = self.k
        self.assertEqual(lhs,rhs,'Multiplication rules fail, ne(ij,k)')
        
    def test_jik(self):
        lhs = self.j*self.i
        rhs = -self.k
        self.assertEqual(lhs,rhs,'Multiplication rules fail, ne(ji,-k)')

    def test_ijji(self):
        lhs = self.i*self.j
        rhs = -(self.j*self.i)
        self.assertEqual(lhs,rhs,'Multiplication rules fail, ne(ij,-ji)')


class NormTestCase(unittest.TestCase):

    def setUp(self):
        
        import quaternion
        
        self.q = quaternion.Quaternion(1,1,1,1)
        
    def test_norm(self):
        lhs = self.q.norm()
        rhs = 2
        self.assertEqual(lhs,rhs,'Norm rule fails')
        
    def test_abs(self):
        lhs = abs(self.q)
        rhs = 2
        self.assertEqual(lhs,rhs,'abs() rule fails')


class DotTestCase(unittest.TestCase):
    
    def setUp(self):
    
        from quaternion import Quaternion, dot
        from random import uniform
        
        lo,hi = -100,100
        a,b,c,d,e,f = (uniform(lo,hi) for i in range(6))
        self.p = Quaternion(0,a,b,c)
        self.q = Quaternion(0,d,e,f)
        self.lhs = dot(self.p,self.q)
        
    def test_dot_id1(self):
        rhs = (self.p.conj() * self.q + self.q.conj() * self.p) * 0.5
        self.assertEqual(self.lhs,rhs,'First dot product identity fails')
        
    def test_dot_id2(self):
        rhs = (self.p * self.q.conj() + self.q * self.p.conj()) * 0.5
        self.assertEqual(self.lhs,rhs,'Second dot product identity fails')        


class CrossTestCase(unittest.TestCase):
    
    def setUp(self):
    
        from quaternion import Quaternion, cross
        from random import uniform
        
        lo,hi = -100,100
        a,b,c,d,e,f = (uniform(lo,hi) for i in range(6))
        self.p = Quaternion(0,a,b,c)
        self.q = Quaternion(0,d,e,f)
        self.lhs = cross(self.p,self.q)
        
    def test_cross(self):
        rhs = (self.p * self.q - self.q.conj() * self.p.conj()) * 0.5
        self.assertEqual(self.lhs,rhs,'Cross product identity fails')


class ProductTestCase(unittest.TestCase):
    
    def setUp(self):
    
        from quaternion import Quaternion, cross, dot
        from random import uniform
        
        self.dot = dot
        self.cross = cross
        
        lo,hi = -100,100
        a,b,c,d,e,f = (uniform(lo,hi) for i in range(6))
        
        self.p = Quaternion(0,a,b,c)
        self.q = Quaternion(0,d,e,f)
        self.lhs = self.p * self.q
        
    def test_product_identity(self):
        s,v = self.p.scalar(), self.p.vector()
        t,w = self.q.scalar(), self.q.vector()
        rhs = s*t - self.dot(v,w) + s*w + t*v + self.cross(v,w)
        self.assertEqual(self.lhs,rhs,'Product identity fails')
        
class NormaliseTestCase(unittest.TestCase):

    def setUp(self):
        
        import quaternion
        
        self.q = quaternion.Quaternion(1,1,1,1)
        self.n = quaternion.Quaternion(0.5,0.5,0.5,0.5)
    
    def test_normalise(self):
        lhs = self.q.normalise()
        rhs = self.n
        self.assertEqual(lhs,rhs,'Normalise fails')
        
class LogExpTestCase(unittest.TestCase):

    def setUp(self):
        
        from quaternion import Quaternion, exp, log
        from random import uniform
        
        self.exp = exp
        self.log = log
        
        lo,hi = -1,1
        a,b,c,d,e,f,g,h = (uniform(lo,hi) for i in range(8))
        
        self.p = Quaternion(a,b,c,d)
        self.q = Quaternion(e,f,g,h)
    
    def test_exp_log(self):
        lhs = self.exp(self.log(self.p))
        rhs = self.p
        self.assertEqual(lhs,rhs,'exp() does not undo log()')
        
    def test_log_exp(self):
        lhs = self.log(self.exp(self.p))
        rhs = self.p
        self.assertEqual(lhs,rhs,'log() does not undo exp()')
