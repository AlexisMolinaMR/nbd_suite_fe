import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
import os
import numpy as np
import streamlit as st
import py3Dmol
from stmol import showmol
import time

def format_converter(uploaded_file):
    '''
    input:  uploaded ligand file

    output: file; SMILES  
    '''
    file_name = uploaded_file.name

    if file_name.endswith('mol'):
        with open('input_file_{}.mol'.format(int(time.time())), 'w') as f:
            f.write(uploaded_file.getvalue().decode("utf-8"))
        
        os.system("obabel input_file_{}.mol -O out_file.smi".format(int(time.time())))

    elif file_name.endswith('sdf'):
        with open('input_file_{}.sdf'.format(int(time.time())), 'w') as f:
            f.write(uploaded_file.getvalue().decode("utf-8"))
              
        os.system("obabel input_file_{}.sdf -O out_file.smi".format(int(time.time())))

                
    elif file_name.endswith('pdb'):
        with open('input_file_{}.pdb'.format(int(time.time())), 'w') as f:
            f.write(uploaded_file.getvalue().decode("utf-8"))
               
        os.system("obabel input_file_{}.pdb -O out_file.smi".format(int(time.time())))

    

def render_mol(xyz, style, spin=False, color='white', color_prot='#3571C3'):
    '''
    input:  pdb file
            style : []
            spin : Bool
            color : str
            color_prot : str;hex

    output: None; Displayer
    '''

    xyzview = py3Dmol.view(width=600,height=600)
    xyzview.addModel(xyz,'pdb')
    xyzview.setStyle({style:{'color':color_prot}})
    xyzview.setBackgroundColor(color)#('0xeeeeee')
    xyzview.zoomTo()

    if spin:
        xyzview.spin(True)
    else:
        xyzview.spin(False)
    showmol(xyzview, height = 500,width=800)

def makeblock(smi):
    '''
    input:  str; SMILES

    output: RDKit molecular block
    '''

    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mblock = Chem.MolToMolBlock(mol)
    return mblock

def render_lig(xyz):
    '''
    input:  pdb file
            style : []
            spin : Bool
            color : str
            color_prot : str;hex

    output: None; Displayer
    '''
    xyzview = py3Dmol.view()#(width=400,height=400)
    xyzview.addModel(xyz,'mol')
    xyzview.setStyle({'stick':{}})
    xyzview.setBackgroundColor('white')
    xyzview.zoomTo()
    showmol(xyzview,height=500,width=500)

def main():

    html_temp = """
    <center>
    <img src="https://www.asebio.com/sites/default/files/2020-01/NBD_highres.png" alt="nbd_logo" style="width:300px;height:120px;">
    </center>
    <div style ="background-color:#a5d0c5">
    <h1 style ="color:white;text-align:center;">PELE Induced Fit Simulation</h1>
    </div>
    <div style ="background-color:white;padding:40px">
    </div>
    """

    st.markdown(html_temp, unsafe_allow_html = True)
           

    st.subheader('Your system')
    system_file = st.file_uploader("Choose a file in PDB format")

    if system_file is not None:
        st.sidebar.subheader('System visualization settings')
        bcolor = st.sidebar.color_picker('Pick a background color', '#ffffff')
        color_prot = st.sidebar.color_picker('Pick a color for your system', '#3571C3',  key=3)
        style = st.sidebar.selectbox('Choose protein display style',['cartoon', 'line','cross','stick','sphere','clicksphere'])
        spin = st.sidebar.checkbox('Make it spin!', value = False)
                
        xyz = system_file.getvalue().decode("utf-8")
        render_mol(xyz, style=style, spin=spin, color=bcolor, color_prot=color_prot)

    st.subheader('Your ligand')
    ligand_file = st.file_uploader("Choose a file in PDB, sdf or SMILES format", key=1)

    ligand_smiles = st.text_input('Or paste your ligand SMILES', value="")
            
    if ligand_file is not None:
        format_converter(ligand_file)
        with open('out_file.smi', 'r') as f:
            ligand = f.read()
        st.sidebar.subheader('Ligand visualization settings')
        bcolor = st.sidebar.color_picker('Pick a background color', '#ffffff')
        blk=makeblock(ligand)
        render_lig(blk)

    elif ligand_smiles != "":

        st.sidebar.subheader('Ligand visualization settings')
        bcolor = st.sidebar.color_picker('Pick a background color', '#ffffff', key=4)
        blk=makeblock(ligand_smiles)
        render_lig(blk)
            
    st.header('Initial parameters')
    st.subheader('Starting coordinates')


    x = st.number_input('Coordinate X')
    y = st.number_input('Coordinate Y')
    z = st.number_input('Coordinate Z')

    radius = st.slider(
            'Choose a selection radius',
            5.0, 20.0, 7.0, 0.5)

    if x and y and z:
        simulation = st.button('Simulate')
            
        if simulation:
            with st.spinner('Wait for it'):
                time.sleep(5)
            st.success('Once the simulation is finished, the results will be sent to the email adress provided.')
        
if __name__=='__main__':
    main()
