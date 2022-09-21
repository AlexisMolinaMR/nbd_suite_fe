import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from scopy.ScoDruglikeness import molproperty
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import streamlit as st
from rdkit.Chem import Draw
import streamlit.components.v1 as components
import shutil
from os.path import exists
import py3Dmol
from stmol import showmol
from rdkit import Chem
from rdkit.Chem import Descriptors
import scopy
from scopy.ScoDruglikeness import molproperty
import matplotlib.pyplot as plt
from rdkit.Chem import Draw
import streamlit.components.v1 as components
import shutil
from os.path import exists
import py3Dmol
from stmol import showmol


def RuleRadar(mol, name,
           prop_kws={'MW':(100,600),'logP':(-3,6),
                     'nHA':(0,12),'nHD':(0,7),'TPSA':(0,180),
                     'nRot':(0,11),'nRing':(0,6),'nC':(3,35),'nHet':(1,15),
                     'HetRatio':(0.1,1.1), 'fChar':(-4,4),'nRing':(0,30),
                     'QEDmean': (0,1), 'SAscore':(1,10)}
               ):
    """ A radar plot positionning compound's values within the selected filter ranges (pale blue and red).
    By default, the `drug-like soft`_ filter ranges are visualized.

    :param mol: molecular
    :type mol: rdkit.Chem.rdchem.Mol
    :return: A radar plot positionning compound's values
    :rtype: matplotlib.figure.Figure

    .. _drug-like soft:
        http://fafdrugs4.mti.univ-paris-diderot.fr/filters.html

    """
    def _disposedNone(dtype,num_a,num_b):
        if dtype == 'MIN':
            NUM = min([num_a,num_b])
            return NUM - 1.2*abs(NUM)
        else:
            NUM = max([num_a,num_b])
            return NUM + 1.2*abs(NUM)

    items = list(prop_kws.keys())
    props = molproperty.GetProperties(mol,items=items)
    num_prop = len(props)

    for item in items:
        assert (prop_kws[item][0] is not None or
                prop_kws[item][1] is not None
                ), "You need to enter at least an upper or lower limit"

    rule_ceil = np.array([prop_kws[item][1]
    if prop_kws[item][1] is not None else _disposedNone('MAX',prop_kws[item][0],props[item])
    for item in items])

    rule_floor = np.array([prop_kws[item][0]
    if prop_kws[item][0] is not None else _disposedNone('MIN',prop_kws[item][1],props[item])
    for item in items])

    props = np.array(list(props.values()))

    bench_floor = np.vstack((props, rule_floor)).min(axis=0)
    bench_floor -= 0.2*bench_floor
    bench_ceil = np.vstack((props, rule_ceil)).max(axis=0)*1.2

    #max-min standarize
    props = (props-bench_floor)/(bench_ceil-bench_floor)
    floor = (rule_floor-bench_floor)/(bench_ceil-bench_floor)
    ceil = (rule_ceil-bench_floor)/(bench_ceil-bench_floor)

    theta = np.linspace(0, 360, num_prop, endpoint=False)
    X_ticks = np.radians(theta)#angle to radian
    X_ticks = np.append(X_ticks,X_ticks[0])
    Y = np.vstack((props,floor,ceil))
    Y = np.hstack((Y, Y[:,0].reshape(3,1)))

    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(10,10))
    ax.plot(X_ticks, Y[0])
    ax.plot(X_ticks, Y[1],color='#339999')
    ax.fill(X_ticks, Y[1], alpha=0.25,color='#339999')
    ax.plot(X_ticks, Y[2],color='#339999')
    ax.fill(X_ticks, Y[2], alpha=0.20,color='#339999')
    ax.set_xticks(X_ticks)
    items.append(items[0])
    ax.set_xticklabels(items, fontsize=10)
    ax.set_yticks([])

    ax.spines['polar'].set_visible(False)
    ax.grid(axis='y')
    ax.set_ylim([0,1])
    for i in [0,0.2,0.4,0.6,0.8,1.0]:
        ax.plot(X_ticks,[i]*(num_prop+1),'-', color='black',lw=0.5)
    ax.set_theta_zero_location('N')
    plt.title('Property profile for {}'.format(name), fontsize=13)
    #plt.tight_layout()
    #plt.show()
    return fig

def makeblock(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mblock = Chem.MolToMolBlock(mol)
    return mblock

def render_mol(xyz):
    xyzview = py3Dmol.view()#(width=400,height=400)
    xyzview.addModel(xyz,'mol')
    xyzview.setStyle({'stick':{}})
    xyzview.setBackgroundColor('white')
    xyzview.zoomTo()
    showmol(xyzview,height=500,width=500)

def compute_descriptors(input_file):

    #in_df = pd.read_csv(input_file)
    desc_df = pd.DataFrame(columns=['canon_smiles', 'inchikey', 'mw', 'hba', 'hbd', 'rtb', 'heteroatoms',
    'aromatic_rings', 'heavy_atoms', 'psa', 'logp', 'logd', 'logsw', 'rgb', 'C',
    'n_stereocenters', 'Fsp3'])

    for i in ['CCCCCc1cc(O)c2c(c1C(=O)O)OC(C)(C)[C@@H]1CCC(C)=C[C@@H]21']:
        l_desc = []

        mol = Chem.MolFromSmiles(i)
        inchikey = Chem.MolToInchi(mol)
        logp = molproperty.CalculateLogP(mol)
        logd = molproperty.CalculateLogD(mol)
        logsw = molproperty.CalculateLogSw(mol)
        rgb=molproperty.CalculateNumRigidBonds(mol)
        C = molproperty.CalculateNumCarbon(mol)
        nster = molproperty.CalculateNumStereocenters(mol)
        fsp3=molproperty.CalculateFsp3(mol)
        hba = Descriptors.NumHAcceptors(mol)
        hbd = Descriptors.NumHDonors(mol)
        rtb = Descriptors.NumRotatableBonds(mol)
        het = Descriptors.NumHeteroatoms(mol)
        aromatic_R = Descriptors.NumAromaticRings(mol)
        hv_at = Descriptors.HeavyAtomCount(mol)
        psa = Descriptors.TPSA(mol)
        mw = Descriptors.ExactMolWt(mol)

        l_desc = [i, inchikey, mw, hba, hbd, rtb, het, aromatic_R, hv_at, psa, logp,
        logd, logsw, rgb, C, nster, fsp3]

        desc_df.loc[len(desc_df)] = l_desc


    return desc_df

@st.cache
def convert_df(df):
    return df.to_csv().encode('utf-8')

def main():

    html_temp = """
    <style>
        #MainMenu {visibility: hidden;}
        input[type="text"], input[type="number"]{border-radius: 0;padding: 0.375rem 0.75rem;font-size: 1rem;font-weight: 400;line-height: 1.5;color: #495057;background-color: #fff;background-clip: padding-box;border: 1px solid #b9b9b8;}
        .css-demzbm {background-color: #006c7b;}
        .css-demzbm:focus {box-shadow: #006c7b80 0px 0px 0px 0.2rem;}
        .st-ds {border-bottom-color: #006c7b;}
        .st-dr {border-top-color: #006c7b;}
        .st-dq {border-right-color: #006c7b;}
        .st-dp {border-left-color: #006c7b;}
        .step-down, .step-up {background-color: transparent;border: 1px solid #b9b9b8;border-radius: 0!important;border-left: 0;}
        .css-1cndplc {height: 38px;}
        .step-down:hover:enabled, .step-down:focus:enabled {background-color: #006c7b;}
        .step-up:hover:enabled, .step-up:focus:enabled {background-color: #006c7b;}
        .css-1cpxqw2 {border-radius: 50px;border-color: #006c7b;color: #006c7b;}
        .css-1cpxqw2:hover {background-color: #006c7b;color: #fff;border-color: #006c7b;}
        .css-1cpxqw2:focus:not(:active) {border-color: #006c7b;color: #006c7b;}
        .css-1cpxqw2:focus {box-shadow: rgb(0 108 123 / 50%) 0px 0px 0px 0.2rem;outline: none;}
        .css-1cpxqw2:hover:focus{color: #fff;}
        .css-1cpxqw2:active {color: rgb(255, 255, 255);border-color: #006c7b;background-color: #006c7b;}
        </style>
    <center>
    <img src="https://www.asebio.com/sites/default/files/2020-01/NBD_highres.png" alt="nbd_logo" style="width:300px;height:120px;">
    </center>
    <div style ="background-color:#a5d0c5">
    <h1 style ="color:white;text-align:center;">NBD Compound Profiling</h1>
    </div>
    <div style ="background-color:white;padding:40px">
    </div>
    """



    st.markdown(html_temp, unsafe_allow_html = True)

    name = st.text_input('Compound name', 'compound')
    smile = st.text_input('Input SMILES', 'CCCCCc1cc(O)c2c(c1C(=O)O)OC(C)(C)[C@@H]1CCC(C)=C[C@@H]21')
    mol = Chem.MolFromSmiles(smile)
    fig_2d = Draw.MolToImage(mol)

    uploaded_file = st.file_uploader("Or choose a file")


    if exists('pdbtemp'):
        shutil.rmtree('pdbtemp')

    if uploaded_file is not None:
        pass

    else:

        if mol is not None:

            if st.button("Generate profile"):
                with st.spinner('Generating...'):


                    fig = RuleRadar(mol,name, {'MW':(100,600),'logP':(-3,6),
                                 'nHA':(0,12),'nHD':(0,7),'TPSA':(0,180),
                                 'nRot':(0,11),'nRing':(0,6),'MaxRing':(0,18),
                                 'nC':(3,35),'nHet':(1,15),'HetRatio':(0.1,1.1),
                                 'fChar':(-4,4), 'QEDmean': (0,1), 'SAscore':(1,10)},
                                 )
                st.success('Done!')
                col1, col2 = st.columns(2)

                with col1:
                    with st.container():
                        st.header("2D structure")
                        st.image(fig_2d)

                with col2:
                    with st.container():
                        st.header("3D structure")
                        blk=makeblock(smile)
                        render_mol(blk)

                prop_df = compute_descriptors(smile)

                st.dataframe(prop_df)


                col3, col4 = st.columns(2)
                
                with col3:
                    st.header("ChEMBL")
                    st.pyplot(fig, True)
                with col4:
                    st.header("ENAMINE")
                    st.pyplot(fig, True)
                
                st.header("ChemistriX")
                st.pyplot(fig, True)



        else:
            st.error('Unvalid molecule')


if __name__=='__main__':
    main()
