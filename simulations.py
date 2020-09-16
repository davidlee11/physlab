import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation,rc
from IPython.display import HTML
from scipy.integrate import odeint


class projectile():
    
    def __init__(self, h, v, s=0, theta=0, dt=0.2):
        self.s, self.h = s, h
        self.v, self.theta = v, theta
        
        self.vx, self.vy = v*np.cos(np.pi*theta/180), v*np.sin(np.pi*theta/180)
        self.dt=dt
        self.time=0
        
        self.tmax=2*self.vy/9.8-self.vy/9.8+np.sqrt((self.vy/9.8)**2+2*self.h/9.8)
        self.frames=int(self.tmax/self.dt)+1
        self.interval=self.dt*1e3
        
        self.fig, self.ax = plt.subplots()
        plt.close()
        
        self.ax.set_xlabel(r'$s (m)$',fontsize=12)
        self.ax.set_ylabel(r'$h (m)$',fontsize=12,rotation=0)
        self.ax.yaxis.set_label_coords(-0.12,0.5)

        self.state = [[], []]
        self.line, = self.ax.plot([],[],'o',color='b')
        self.text = self.ax.text(0.01,1.01,r'$time : {:.2f}s,\quad s : {:.2f}m$'.format(self.time,self.s)
                                 ,fontsize=12,transform=self.ax.transAxes)

    def init_func(self):
        self.line.set_data([],[])
        return (self.line,)

    def animate(self,i):
        if self.h>=0:
            self.state[0].append(self.s)
            self.state[1].append(self.h)

        self.line.set_data(self.state[0],self.state[1])
        self.ax.set_xlim(-1,1+self.state[0][-1]*1.1)
        self.ax.set_ylim(-1,max(self.state[1])*1.1+1)
        self.text.set_text(r'$time : {:.2f}s,\quad s : {:.2f}m$'.format(self.time,self.s))
        
        self.s += self.vx*self.dt
        self.h += self.vy*self.dt-0.5*9.8*self.dt**2
        self.vy -= 9.8*self.dt
        
        self.time += self.dt

        return (self.line,)
    
    def call(self):
        anim=animation.FuncAnimation(self.fig,
                                     self.animate,
                                     init_func=self.init_func,
                                     frames=self.frames,
                                     interval=self.interval,
                                     blit=True)
        
        rc('animation',html='jshtml')
        
        return anim
    

class momentum_conservation():
    
    def __init__(self, m1=2, v1=3, m2=1, v2=1, dt=0.1):
        self.x1, self.y1, self.m1, self.v1 = 0, 0.5, m1, v1
        self.x2, self.y2, self.m2, self.v2 = 10, 0.5, m2, v2
        
        self.momentum=self.m1*self.v1+self.m2*self.v2
        self.dt=dt
        self.time=0
        
        self.tmax=2*10/abs(self.v1-self.v2)
        self.frames=int(self.tmax/self.dt)+1
        self.interval=self.dt*1e3
        
        self.fig, self.ax = plt.subplots()
        plt.close()
        
        self.state1 = [[], [], [], []]
        self.state2 = [[], [], [], []]
        self.obj1, = self.ax.plot([],[],'o',ms=10,color='b')
        self.obj2, = self.ax.plot([],[],'o',ms=10,color='r')
        self.arrow = [self.ax.arrow(0,0,0,0),self.ax.arrow(0,0,0,0)]
        self.text = self.ax.text(0.01,1.01,r'$p : {:.2f}kg m/s,\quad v_1 : {:.2f}m/s,\quad v_2 : {:.2f}m/s$'.format(self.momentum,self.v1,self.v2)
                                 ,fontsize=12,transform=self.ax.transAxes)

    def init_func(self):
        self.obj1.set_data([],[])
        self.obj2.set_data([],[])
        return (self.obj1, self.obj2,)

    def animate(self,i):
                       
        if abs(self.x1-self.x2)<=1:
            v1_f=(self.m1-self.m2)/(self.m1+self.m2)*self.v1+2*self.m2/(self.m1+self.m2)*self.v2
            v2_f=2*self.m1/(self.m1+self.m2)*self.v1-(self.m1-self.m2)/(self.m1+self.m2)*self.v2
            self.v1,self.v2=v1_f,v2_f
            
            
        self.obj1.set_data(self.x1,self.y1)
        self.obj2.set_data(self.x2,self.y2)
        self.arrow[0].remove()
        self.arrow[1].remove()
        self.arrow[0]=self.ax.arrow(self.x1,self.y1,self.v1,0,
                                    head_width=0.2,head_length=0.5,color='b',lw=2,zorder=1)
        self.arrow[1]=self.ax.arrow(self.x2,self.y2,self.v2,0,
                                    head_width=0.2,head_length=0.5,color='r',lw=2,zorder=1)
        
        
        self.ax.set_xlim(-1,30)
        self.ax.set_ylim(-1,3)
        self.text.set_text(r'$p : {:.2f}kg \cdot m/s,\quad v_1 : {:.2f}m/s,\quad v_2 : {:.2f}m/s$'.format(self.momentum,self.v1,self.v2))
        
        self.x1 += self.v1*self.dt
        self.x2 += self.v2*self.dt
        
        self.time += self.dt

        return (self.obj1, self.obj2,)
    
    def call(self):
        anim=animation.FuncAnimation(self.fig,
                                     self.animate,
                                     init_func=self.init_func,
                                     frames=self.frames,
                                     interval=self.interval,
                                     blit=True)
        
        rc('animation',html='jshtml')
        
        return anim
    

class energy_conservation():
    
    def __init__(self, h, v=0, m=1, dt=0.2):
        self.h = h
        self.v = v
        self.m = m
        
        self.dt=dt
        self.time=0
        
        self.tmax=np.sqrt(2*self.h/9.8)
        self.frames=int(self.tmax/self.dt)+1
        self.interval=self.dt*1e3
        
        self.fig, self.ax = plt.subplots()
        plt.close()
        
        self.ax.set_ylabel(r'$h (m)$',fontsize=12,rotation=0)
        self.ax.yaxis.set_label_coords(-0.12,0.5)

        self.state = [[], [], [], []]
        self.line, = self.ax.plot([],[],'o',color='k')
        self.bar=[self.ax.bar(0,0,color='r'),self.ax.bar(0,0,color='g'),self.ax.bar(0,0,color='b')]
        self.text = self.ax.text(0.01,1.01,r'$h : {:.2f}m,\quad v : {:.2f}m/s$'.format(self.h,self.v)
                                 ,fontsize=12,transform=self.ax.transAxes)
        
        self.ax.set_xticklabels(('', '', 'Kinetic', 'Potentional', 'Mechanical'))

    def init_func(self):
        self.line.set_data([],[])
        return (self.line,)

    def animate(self,i):
        T=0.5*self.m*self.v**2
        V=self.m*9.8*self.h
        self.state[0].append(self.h)
        self.state[1].append(T*0.1)
        self.state[2].append(V*0.1)
        self.state[3].append(T*0.1+V*0.1)
        
        
        self.line.set_data(0,self.state[0])
        self.bar[0].remove()
        self.bar[1].remove()
        self.bar[2].remove()
        self.bar[0]=self.ax.bar(1,self.state[1][-1],color='r')
        self.bar[1]=self.ax.bar(2,self.state[2][-1],color='g')
        self.bar[2]=self.ax.bar(3,self.state[3][-1],color='b')
        self.ax.set_xlim(-1,4)
        self.ax.set_ylim(0,max(self.state[0])*1.1+1)
        self.text.set_text(r'$h : {:.2f}m,\quad v : {:.2f}m/s,\quad E : {:.2f}J$'.format(self.h,self.v,T+V))
        
        self.h += self.v*self.dt-0.5*9.8*self.dt**2
        self.v -= 9.8*self.dt
        
        self.time += self.dt

        return (self.line,)
    
    def call(self):
        anim=animation.FuncAnimation(self.fig,
                                     self.animate,
                                     init_func=self.init_func,
                                     frames=self.frames,
                                     interval=self.interval,
                                     blit=True)
        
        rc('animation',html='jshtml')
        
        return anim
    

class ballistic_pendulum():
    
    def __init__(self, m=30e-3, M=870e-3, v=100, l=2, dt=0.1):
        self.x1, self.y1, self.m1, self.v1 = -5, -l, m, v
        self.x2, self.y2, self.m2, self.v2 = 0, -l, M, 0
        self.l = l
        
        self.V=0
        self.omega=self.V/self.l
        self.theta=0
        
        self.dt=dt
        self.time=0
        
        self.tmax=20
        self.frames=int(self.tmax/self.dt)+1
        self.interval=self.dt*1e3
        
        self.cnt=0
        
        self.fig, self.ax = plt.subplots(figsize=(12,4))
        plt.close()
        
        self.state1 = [[], [], [], []]
        self.state2 = [[], [], [], []]
        self.pivot = self.ax.plot(0,0,'o',ms=10,color='k')
        self.obj1, = self.ax.plot([],[],'o',ms=6,color='purple',zorder=3)
        self.obj2, = self.ax.plot([],[],'o',ms=20,color='r',zorder=2)
        self.line, = self.ax.plot([],[],lw=2,color='purple',zorder=1)
        self.text = self.ax.text(0.01,1.01,r'$time : {:.2f}s, \theta : {:.2f}^\degree$'.format(self.time,180*self.theta/np.pi)
                                 ,fontsize=12,transform=self.ax.transAxes)

    def init_func(self):
        self.obj1.set_data([],[])
        self.obj2.set_data([],[])
        self.line.set_data([],[])
        return (self.obj1, self.obj2, self.line,)

    def animate(self,i):
                       
        if self.x1==0 and self.cnt==0:
            self.V=self.m1/(self.m1+self.m2)*self.v1
            self.cnt=1
        
        self.obj1.set_data(self.x1,self.y1)
        self.obj2.set_data(self.x2,self.y2)
        self.line.set_data([0,self.x2],[0,self.y2])
                
        self.ax.set_xlim(-7,5)
        self.ax.set_ylim(-3,1)
        self.text.set_text(r'$time : {:2.2f}s, \theta : {:.2f}^\degree$'.format(self.time,180*self.theta/np.pi))
               
        self.omega=self.V/self.l
        
        self.theta, self.omega = odeint(self.pend, [self.theta,self.omega], [0,self.dt])[-1]
        self.V=self.l*self.omega
        
        self.x1=self.x2=self.l*np.sin(self.theta)
        self.y1=self.y2=-self.l*np.cos(self.theta)
               
        self.time += self.dt

        return (self.obj1, self.obj2, self.line,)
    
    def pend(self,y,t):
        theta, omega = y
        dydt = [omega, -(9.8/self.l)*np.sin(theta)]
        return dydt
    
    def call(self):
        anim=animation.FuncAnimation(self.fig,
                                     self.animate,
                                     init_func=self.init_func,
                                     frames=self.frames,
                                     interval=self.interval,
                                     blit=True)
        
        rc('animation',html='jshtml')
        
        return anim


class wave_properties():
    
    def __init__(self, k, w, phi=0, dt=0.1):
        self.k, self.w, self.phi = k, w, phi*np.pi/180

        self.dt=dt
        self.time=0
        
        self.tmax=2*np.pi/self.w
        self.frames=int(self.tmax/self.dt)+1
        self.interval=self.dt*1e3
        
        self.fig, self.ax = plt.subplots(figsize=(15,4))
        plt.close()
        
        self.ax.set_xlabel(r'$x (m)$',fontsize=12)
        self.ax.set_ylabel(r'$y (m)$',fontsize=12,rotation=0)
        self.ax.yaxis.set_label_coords(-0.12,0.5)

        self.line, = self.ax.plot([],[],'o',color='k')
        self.point, = self.ax.plot([],[],'o',color='r')
        self.pline, = self.ax.plot([],[],color='r')

        self.cpoint, = self.ax.plot([],[],'o',color='r')
        self.cpline, = self.ax.plot([],[],color='k')
        self.projline, = self.ax.plot([],[],color='r')
        self.circle = self.ax.plot(np.cos(np.linspace(0,2*np.pi,51))-2.5,
                                   np.sin(np.linspace(0,2*np.pi,51)),':',color='k')
        
        self.ltheta, = self.ax.plot([],[],color='b')
        self.ttheta = self.ax.text(-2.5,1.1,r'$\theta : {:.4g}\degree$'.format(self.phi*180/np.pi),
                                   fontsize=12)
        self.linebtw, = self.ax.plot([],[],color='b')
        
        self.ax.axhline(0,lw=1,color='k')
        self.ax.axvline(0,lw=1,color='k')
        self.ax.text(np.pi/self.k,1.5,r'$\lambda : {:.4g}m$'.format(2*np.pi/self.k),fontsize=12,ha='center',va='bottom')
        self.ax.text(0.01,0.95,r'$f : {:.4g}Hz$'.format(self.w/(2*np.pi))+'\n'+r'$T : {:.4g}s$'.format(self.tmax),
                     fontsize=12,va='top',transform=self.ax.transAxes)
        self.ax.annotate('',(2*np.pi/self.k,1.5),(0,1.5),arrowprops=dict(arrowstyle='<->'),color='g')
        self.ax.plot([-2.5,np.cos(self.phi)-2.5],[0,np.sin(self.phi)],lw=1,color='k')

        self.text = self.ax.text(0.01,1.01,r'$time : {:.2f}s'.format(self.time),
                                 fontsize=12,transform=self.ax.transAxes)

    def init_func(self):
        self.line.set_data([],[])
        self.point.set_data([],[])
        self.pline.set_data([],[])
        return (self.line, self.point, self.pline)

    def animate(self,i):
        x=np.linspace(0,10,int(51*self.k))
        y=np.sin(self.k*x-self.w*self.time+self.phi)
        cx=np.cos(-self.w*self.time+self.phi)-2.5
        cy=np.sin(-self.w*self.time+self.phi)

        self.line.set_data(x,y)
        self.point.set_data(x[0],y[0])
        self.pline.set_data([x[0],x[0]],[0,y[0]])
        self.cpoint.set_data(cx,cy)
        self.cpline.set_data([-2.5,cx],[0,cy])
        self.projline.set_data([cx,cx],[0,cy])
        self.ltheta.set_data(0.25*np.cos(-np.arange(0,self.w*self.time,0.1)+self.phi)-2.5,
                             0.25*np.sin(-np.arange(0,self.w*self.time,0.1)+self.phi))
        self.ttheta.set_text(r'$\theta : {:.4g}\degree$'.format(((self.w*self.time+self.phi)*180/np.pi)%360))
        self.linebtw.set_data([cx,x[0]],[cy,y[0]])

        self.ax.imshow(np.array([y for f in range(2)]),cmap='gray',extent=[0,10,-1.75,-1.25])

        self.ax.set_xlim(-5,10)
        self.ax.set_ylim(-2,2)
        self.text.set_text(r'$time : {:.2f}s$'.format(self.time))
        
        self.time += self.dt

        return (self.line, self.point, self.pline)
    
    def call(self):
        anim=animation.FuncAnimation(self.fig,
                                     self.animate,
                                     init_func=self.init_func,
                                     frames=self.frames,
                                     interval=self.interval,
                                     blit=True)
        
        rc('animation',html='jshtml')
        
        return anim
