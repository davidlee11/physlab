import numpy as np
import plotly.offline as py
import plotly.graph_objects as go
from ipywidgets import interact,IntSlider,FloatSlider

class diffraction():

    def single(self):
        style={'description_width':'initial'}
        @interact(wave_length=IntSlider(value=600,min=400,max=700,step=1,description='파장 '+chr(0x03BB)+' (nm)',style=style),
                  slit_width=FloatSlider(value=0.5,min=0.1,max=1,step=0.1,description='슬릿 넓이 w (mm)',style=style),
                  screen_distance=IntSlider(value=100,min=50,max=200,step=1,description='스크린 거리 L (cm)',style=style))

        def update(wave_length,slit_width,screen_distance):
            scale=1e3
            x=np.linspace(-2e-3,2e-3,1000)
            y=self.single_slit_diffraction(x,wave_length,slit_width,screen_distance)
            fig=go.FigureWidget()
            fig.update_layout(autosize=False,width=800,height=400,
                              title='단일 슬릿 시뮬레이션',
                              xaxis=dict(title='Distance (mm)'),
                              yaxis=dict(title='Intensity',tickvals=[0,0.5,1]))
            fig.add_trace(go.Heatmap(z=[y for f in range(2)],x=x*scale,y=[1.25,1.5,1.75],colorscale='Greys',reversescale=True))
            fig.add_scatter(x=x*scale,y=y)
            fig.show()

    def double(self):
        style={'description_width':'initial'}
        @interact(wave_length=IntSlider(value=600,min=400,max=700,step=1,description='파장 '+chr(0x03BB)+' (nm)',style=style),
                  slit_width=FloatSlider(value=0.5,min=0.1,max=1,step=0.1,description='슬릿 넓이 w (mm)',style=style),
                  slit_distance=FloatSlider(value=3,min=1,max=5,step=0.1,description='슬릿 간격 d (mm)',style=style),
                  screen_distance=IntSlider(value=100,min=50,max=200,step=1,description='스크린 거리 L (cm)',style=style))

        def update(wave_length,slit_width,slit_distance,screen_distance):
            scale=1e3
            x=np.linspace(-2e-3,2e-3,1000)
            y=self.double_slit_diffraction(x,wave_length,slit_width,slit_distance,screen_distance)
            fig=go.FigureWidget()
            fig.update_layout(autosize=False,width=800,height=400,
                              title='이중 슬릿 시뮬레이션',
                              xaxis=dict(title='Distance (mm)'),
                              yaxis=dict(title='Intensity',tickvals=[0,0.5,1]))
            fig.add_trace(go.Heatmap(z=[y for f in range(2)],x=x*scale,y=[1.25,1.5,1.75],colorscale='Greys',reversescale=True))
            fig.add_scatter(x=x*scale,y=y)
            fig.show()

    def single_slit_diffraction(self,x,wave_length,slit_width,screen_distance):
        l=wave_length*1e-9
        a=slit_width*1e-3
        L=screen_distance*1e-2
        T=l*L/a
        Ix=(np.sin(np.pi*a*x/(l*L))/(np.pi*a*x/(l*L)))**2
        return np.where(abs(x)<T,Ix,0)    

    def double_slit_diffraction(self,x,wave_length,slit_width,slit_distance,screen_distance):
        l=wave_length*1e-9
        a=slit_width*1e-3
        d=slit_distance*1e-3
        L=screen_distance*1e-2
        T=l*L/a
        Ix=(np.sin(np.pi*a*x/(l*L))/(np.pi*a*x/(l*L)))**2* \
        np.cos(np.pi*d*x/(l*L))**2
        return np.where(abs(x)<T,Ix,0)