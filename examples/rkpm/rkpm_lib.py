from popcorn import *
import popcorn.libs.gsl_pop as gsl_pop
from popcorn.functional import *
#
# The basis and influence function
#
Polys = { 1: lambda x,y : Matrix([ 1.0, x, y ]),
          2: lambda x,y : Matrix([ 1.0, x, y, x**2, x*y, y**2 ])
        }
Weights = {
    'const' :  lambda r : 1.0,
    'cubic' :  lambda r : (1.0-r)**3,
    'cubicspline':lambda r : Piecewise(
        (2.0/3.0 - 4.0*r**2 + 4.0* r**3, Le(r,0.5)  ),
        (4.0/3.0 - 4.0*r + 4.0*r**2 - 4.0/3.0*r**3, Lt(r,1.0)),
        (0.0,True) ),
    'quarticA':lambda r : (1.0 - 6.0*r**2 + 8.0*r**3 - 3.0*r**4),
    'quarticB':lambda r : (2.0/3.0 - 9.0/2.0*r**2 + 19.0/3.0*r**3 - 5.0/2.0*r**4),
    }
gdim=2
II, JJ, l_edge = symbols('II JJ l_edge')
Npt = l_edge -1
class RKPM_Basis():
    def __init__(self, i_x, y, SupRad, P, w, grad="MLS"):
        self.i_x = i_x

        x0, xI, xJ = self.i_x.Vertex_Handles(0,II,JJ)

        # Geoemetry
        rI = xI - y
        rJ = xJ - y
        rIabs = sqrt((rI.T*rI)[0,0])
        rJabs = sqrt((rJ.T*rJ)[0,0])
        self.y = y
        
        self.P0 = P( 0,0 )
        sx,sy = symbols('x y')
        self.grad_P0 = P( sx,sy ).jacobian([sx,sy]).subs(sx,0).subs(sy,0)

        self.PI = P( *rI )
        self.PJ = P( *rJ )
        self.wI = w( rIabs / SupRad )
        self.wJ = w( rJabs / SupRad )

        # self.wI = w(abs(rI[0])/SupRad)*w(abs(rI[1])/SupRad)
        # self.wJ = w(abs(rJ[0])/SupRad)*w(abs(rJ[1])/SupRad)
        self.NP = len(self.P0)

        M_expr = self.wI * self.PI * self.PI.T
        self.M_expr = M_expr
        
        self.pv_M     = PopcornVariable('rkpm_M', self.NP, 2 )
        self.pv_c   = PopcornVariable('rkpm_c', self.NP, 1 )
        pv_P0  = PopcornVariable('rkpm_P0', self.NP, 1)

        self.init_prgm = [
            self.pv_M,
            Loop(II,0,Npt, [
                Asgn( self.pv_M,   M_expr, '+=')
            ]),
            # DebugPrint(pv_M),
            gsl_pop.mat_lu(self.pv_M),
            # DebugPrint(pv_M),
            pv_P0,
            self.pv_c,
            Asgn(pv_P0, self.P0 ),
            gsl_pop.mat_lu_solve(self.pv_M, self.pv_c, pv_P0),
        ]
        self.close_prgm = [
            gsl_pop.mat_lu_cleanup(self.pv_M)
        ]

        self.NI = self.wI*self.pv_c.as_matrix().dot(self.PI)
        self.NJ = self.wJ*self.pv_c.as_matrix().dot(self.PJ)

        if grad=="MLS":
            self.MLS_Grad()
        else:
            self.Synchronized_Grad()

        
    def Synchronized_Grad(self):
        self.pv_d0_c = PopcornVariable('rkpm_d0_c', self.NP, 1 )
        self.pv_d1_c = PopcornVariable('rkpm_d1_c', self.NP, 1 )

        pv_d0_P0  = PopcornVariable('rkpm_d0_P0', self.NP, 1)
        pv_d1_P0  = PopcornVariable('rkpm_d1_P0', self.NP, 1)

        self.calc_d_c_prgm = [
            self.pv_d0_c,
            self.pv_d1_c,
            pv_d0_P0,
            pv_d1_P0,

            Asgn(pv_d0_P0, self.grad_P0[:,0] ),
            Asgn(pv_d1_P0, self.grad_P0[:,1] ),

            gsl_pop.mat_lu_solve(self.pv_M, self.pv_d0_c, pv_d0_P0),
            gsl_pop.mat_lu_solve(self.pv_M, self.pv_d1_c, pv_d1_P0)
        ]
        
        self.init_prgm += self.calc_d_c_prgm
        
        self.grad_NI = Matrix([
            self.pv_d0_c.as_matrix().T * self.wI*self.PI,
            self.pv_d1_c.as_matrix().T * self.wI*self.PI
        ])
        self.grad_NJ = Matrix([
            self.pv_d0_c.as_matrix().T * self.wJ*self.PJ,
            self.pv_d1_c.as_matrix().T * self.wJ*self.PJ
        ])

        
    def MLS_Grad(self):
        d_M_expr = [ self.M_expr.diff(z) for z in self.y ]
        pv_c_dot_d_M = [ PopcornVariable('rkpm_c_dot_d{0}_M'.format(i),self.NP,1)
                       for i in xrange(gdim) ]
        pv_c_dotdM_Mi = [ PopcornVariable('rkpm_c_dot_d{0}_M_Mi'.format(i),self.NP,1)
                              for i in xrange(gdim) ]
        self.mls_grad_prgm = pv_c_dot_d_M +  pv_c_dotdM_Mi + \
          [
              Loop(II,0,Npt, [
                  Asgn( pv_c_dot_d_M[i], d_M_expr[i] * self.pv_c.as_matrix(), '+=')
                ])
            for i in xrange(gdim) ] +\
          [ gsl_pop.mat_lu_solve(self.pv_M, pv_c_dotdM_Mi[i] ,pv_c_dot_d_M[i])
                for i in xrange(gdim) ]
        self.grad_NI = Matrix([
            - pv_c_dotdM_Mi[i].as_matrix().dot( self.wI * self.PI ) +
            self.pv_c.as_matrix().dot( (self.wI * self.PI).diff( self.y[i] ) )
            for i in xrange(gdim) ])
        self.grad_NJ = Matrix([
            - pv_c_dotdM_Mi[i].as_matrix().dot( self.wJ * self.PJ ) +
            self.pv_c.as_matrix().dot( (self.wJ * self.PJ).diff( self.y[i] ) )
            for i in xrange(gdim) ])

        #from IPython import embed ; embed()
        self.init_prgm += self.mls_grad_prgm

        
    def Interpolate(self, name, i_u, o_uy):
        uI = i_u.Vertex_Handle(II)
        Kernel(name,
               listing = self.init_prgm + [
                   Loop(II,0,Npt,[
                       Asgn(o_uy, self.NI * uI, '+=')
                   ]),
               ] + self.close_prgm
        )
    
    def Integral(self, name, P,tu,u, i_u,ospace, make_K=True):
        o_R = Output('R', [ospace], 1)
        o_K = Output('K', [ospace], 2)
        R = gateaux(P, tu,   Matrix([self.NI]),self.grad_NI)
        #R.simplify()
        K = gateaux(R,  u,   Matrix([self.NJ]),self.grad_NJ)
        uI = i_u.Vertex_Handle(II)
        Kernel(name,
               listing = self.init_prgm + [
                   u,
                   grad(u),
                   Loop(II,0,Npt,[
                       Asgn(u,      self.NI * uI, '+='),
                       Asgn(grad(u), self.grad_NI * uI, '+=')
                   ]),
                   Loop(II,0,Npt,[
                       Asgn(o_R.View((II,)), Matrix([R]),"+="),
                       Loop(JJ,0,Npt,[
                           Asgn(o_K.View((II,JJ)), Matrix([K]),"+=")
                       ])
                   ]),
               ] + self.close_prgm
        )
