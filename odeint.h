class Ode_Int {
public:
  int ignore, kount;
  double dxsav, minstep;
  void init(int n);
  void tidy(void);
  void go(double x1, double x2, double xstep,
	  double eps, void (*derivs)(double, double[],double[]));
  void go_simple(double x1, double x2, int nstep,
		 void (*derivs)(double, double[],double[]));
  void go_scale(double x1, double x2, double step,
		 void (*derivs)(double, double[],double[]));
  void set_bc(int n, double num);
  double get_x(int i);
  double get_y(int n, int i);
  double get_d(int n, int i);
  double *xp, **yp;

private:
  double **dydxp,*hstr,*ystart;
  int kmax,nok,nbad,nvar;
  void rkck(double y[], double dydx[], int n, double x, double h,
	    double yout[],
	    double yerr[], void (*derivs)(double, double[], double[]));
  void rkqs(double y[], double dydx[], int n, double *x, double htry, 
	    double eps,	double yscal[], double *hdid, double *hnext,
	    void (*derivs)(double, double[], double[]));
  void odeint(double ystart[], int nvar, double x1, double x2, double eps, 
	      double h1,double hmin, int *nok, int *nbad,
	      void (*derivs)(double, double [], double []));
#define float double
  void rk4(float y[], float dydx[], int n, float x, float h, float yout[],
	     void (*derivs)(float, float [], float []));
  void rkdumb(float vstart[], int nvar, float x1, float x2, int nstep,
	void (*derivs)(float, float [], float []));
  void rkscale(float vstart[], int nvar, float x1, float x2, float h1,
	void (*derivs)(float, float [], float []));
#undef float  

};
