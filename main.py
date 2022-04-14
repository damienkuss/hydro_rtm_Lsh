import streamlit as st
import numpy as np
import pandas as pd
from PIL import Image
from utils_lsh import hydro_lsh, read_markdown_file, convert_df
from hydro_rtm import h_critique_mu, Hs_critique_mu, h_ferguson_mu




logo = Image.open('river.png')
st.set_page_config(layout="centered", page_icon=logo, page_title="Calcul L sur h")
st.image("desktop.png")



intro_markdown=read_markdown_file("intro.md")
st.markdown(intro_markdown,unsafe_allow_html=True)
st.markdown("---")

input_sup = st.container()
input_med = st.container()
input_inf = st.container()
graph = st.container()
resultats = st.container()
table = st.container()


st.sidebar.markdown("## Données hydrauliques")
Q=st.sidebar.number_input(u'Débit (m3/s)',value=50.,min_value=0.5,step=10.)
loi_frottement = st.sidebar.selectbox('Loi de frottement',('Régime critique','VPE Ferguson (2007)'))
desactiv_param = loi_frottement=='VPE Ferguson (2007)'
if desactiv_param:
    Pente = st.sidebar.number_input('Pente (m/m)',value=0.02,min_value=0.001,step=0.005,format="%.3f")
    D84 = st.sidebar.number_input('D84 (m)',value=0.1,min_value=0.01,step=0.01)
else :
    Pente = 0.02
    D84=0.1


st.sidebar.markdown("## Données Graph")
side1,side2=st.sidebar.columns(2)
binf_calcul=side1.slider("Largeur min", min_value=0.0000001, max_value=100., value=5., step=5.)
bsup_calcul=side2.slider("Largeur max", min_value=5., max_value=200., value=60., step=5.)

lignes_rappel = st.sidebar.checkbox('Lignes de rappel L sur h',value=True)
if lignes_rappel:
    Lshinf = st.sidebar.number_input('L/h inf',value=15.,min_value=1.,step=5.)
    Lshsup = st.sidebar.number_input('L/h sup',value=40.,min_value=1.,step=5.)
else:
    Lshinf = 15.
    Lshsup = 40.



with graph:
    data_trace = st.multiselect(
     'Séries de données à tracer',
     ['hauteur', 'vitesse', 'charge','L sur h'],
     ['hauteur', 'charge','vitesse','L sur h'])

    calc_hydro=hydro_lsh(Q,binf_calc=binf_calcul,bsup_calc=bsup_calcul,step=0.5,Lshinf=Lshinf,Lshsup=Lshsup)
    calc_hydro.set_d84(D84)
    calc_hydro.set_S(Pente)
    calc_hydro.loi_frottement = loi_frottement 
    fig=calc_hydro.make_fig(data_trace,lignes_rappel)
    st.pyplot(fig)


with resultats:
    st.header('Résultats pour valeur de L sur h spécifique')
    col1_2, col2_2 = st.columns(2)
    col0, col1, col2, col3 = st.columns(4)
    Lsh=col1_2.number_input('Indiquez une valeur de ratio L/h',value=15.,min_value=0.5,step=1.,help='Renseignez la valeur et validez avec touche Entrée')
    calc_hydro.set_bsh(Lsh)
    col0.metric("Largeur", f"{round(calc_hydro.b,1)} (m)")
    col1.metric("Hauteur d'eau", f"{round(calc_hydro.h,2)} (m)")
    col2.metric("Charge", f"{round(calc_hydro.H,2)} (m)")
    col3.metric("Vitesse", f"{round(calc_hydro.u,1)} (m/s)")

with table:
    st.header('Tableau de données')
    data=np.transpose([[round(i,2) for i in calc_hydro.b_range],[round(i,1) for i in calc_hydro.bsh_range],
                      [round(i,2) for i in calc_hydro.h_range],
                      [round(i,1) for i in calc_hydro.u_range],
                      [round(i,2) for i in calc_hydro.H_range],
                      np.ones(len(calc_hydro.b_range))*Q])
    pd.options.display.float_format = '${:,.2f}'.format
    df=pd.DataFrame(data,columns=['b','b/h','h','u','H','Q'])
    df
    csv = convert_df(df)

    st.download_button(
    "Télécharger le tableau",
    csv,
    "file_L_sur_h.csv",
    "text/csv",
    key='download-csv'
    )

# Remove whitespace from the top of the page and sidebar
st.markdown("""
        <style>
               .css-18e3th9 {
                    padding-top: 0rem;
                    padding-bottom: 10rem;
                    padding-left: 5rem;
                    padding-right: 5rem;
                }
               .css-1d391kg {
                    padding-top: 3.5rem;
                    padding-right: 1rem;
                    padding-bottom: 3.5rem;
                    padding-left: 1rem;
                }
        </style>
        """, unsafe_allow_html=True)
