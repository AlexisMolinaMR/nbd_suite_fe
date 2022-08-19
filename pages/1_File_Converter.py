import streamlit as st
import os
import time

@st.cache
def convert_df(df):
    return df.to_csv().encode('utf-8')

def color_pred(val):
    color = 'lightgreen' if val >= 0.5 else 'lightred'
    return f'background-color: {color}'

def main():

    html_temp = """
    <center>
    <img src="https://www.asebio.com/sites/default/files/2020-01/NBD_highres.png" alt="nbd_logo" style="width:300px;height:120px;">
    </center>
    <div style ="background-color:#a5d0c5">
    <h1 style ="color:white;text-align:center;">File format converter</h1>
    </div>
    <div style ="background-color:white;padding:40px">
    </div>
    """


    sidebar_html = '''
    <h6 style ="color:Grey;text-align:center;">Supported conversions:</h6>
    <ul style="list-style-type:circle;">
    <li>SMILES 2 MOL</li>
    <li>SMILES 2 SDF</li>
    <li>SMILES 2 PDB</li>
    </ul>
    <p><span class="small">v0.1. - 032022</span></p>
    '''

    st.markdown(html_temp, unsafe_allow_html = True)
    add_selectbox = st.sidebar.markdown(sidebar_html, unsafe_allow_html = True)

    uploaded_file = st.file_uploader("Choose a file")
    print(uploaded_file)
    file_in = st.text_input('Or paste input file')
    
    if file_in != "":
        uploaded_file = file_in
    
    

    if uploaded_file:
        file_name = uploaded_file.name

        st.header('Molecule conversions')
        st.subheader("Choose output format")

       # if st.button("Convert"):
        if file_name.endswith('smi'):
            mol_out = st.checkbox('MOL', key=6)
            sdf_out = st.checkbox('SDF', key=7)
            pdb_out = st.checkbox('PDB', key=8)

            if st.button("Convert"):
                with open('input_file_{}.smi'.format(int(time.time())), 'w') as f:
                    f.write(uploaded_file.getvalue().decode("utf-8"))

                if mol_out:
                    os.system("obabel input_file_{}.smi -O out_file.mol --gen3D".format(int(time.time())))

                    with open('out_file.mol', 'r') as f:
                        st.download_button(label="Download file as mol", data=f, file_name='out_file.mol')

                if sdf_out:
                    os.system("obabel input_file_{}.smi -O out_file.sdf --gen3D".format(int(time.time())))

                    with open('out_file.sdf', 'r') as f:
                        st.download_button(label="Download file as SDF", data=f, file_name='out_file.sdf')

                if pdb_out:
                    os.system("obabel input_file_{}.smi -O out_file.pdb --gen3D".format(int(time.time())))

                    with open('out_file.pdb', 'r') as f:
                        st.download_button(label="Download file as PDB", data=f, file_name='out_file.pdb')

        elif file_name.endswith('mol'):
            smiles_out = st.checkbox('SMILES', key=5)
            sdf_out = st.checkbox('SDF', key=27)
            pdb_out = st.checkbox('PDB', key=28)
            
            if st.button("Convert"):

                with open('input_file_{}.mol'.format(int(time.time())), 'w') as f:
                    f.write(uploaded_file.getvalue().decode("utf-8"))
                if smiles_out:
                    os.system("obabel input_file_{}.mol -O out_file.smi".format(int(time.time())))

                    with open('out_file.smi', 'r') as f:
                        st.download_button(label="Download file as SMILES", data=f, file_name='out_file.smi')
                
                if sdf_out:
                    os.system("obabel input_file_{}.mol -O out_file.sdf".format(int(time.time())))

                    with open('out_file.sdf', 'r') as f:
                        st.download_button(label="Download file as SDF", data=f, file_name='out_file.sdf')

                if pdb_out:
                    os.system("obabel input_file_{}.mol -O out_file.pdb".format(int(time.time())))

                    with open('out_file.pdb', 'r') as f:
                        st.download_button(label="Download file as PDB", data=f, file_name='out_file.pdb')

        elif file_name.endswith('sdf'):
            smiles_out = st.checkbox('SMILES', key=35)
            mol_out = st.checkbox('MOL', key=36)
            pdb_out = st.checkbox('PDB', key=38)
            if st.button("Convert"):
                with open('input_file_{}.sdf'.format(int(time.time())), 'w') as f:
                    f.write(uploaded_file.getvalue().decode("utf-8"))
                if smiles_out:
                    os.system("obabel input_file_{}.sdf -O out_file.smi".format(int(time.time())))

                    with open('out_file.smi', 'r') as f:
                        st.download_button(label="Download file as SMILES", data=f, file_name='out_file.smi')
                
                if mol_out:
                    os.system("obabel input_file_{}.sdf -O out_file.mol".format(int(time.time())))

                    with open('out_file.mol', 'r') as f:
                        st.download_button(label="Download file as SDF", data=f, file_name='out_file.mol')

                if pdb_out:
                    os.system("obabel input_file_{}.sdf -O out_file.pdb".format(int(time.time())))

                    with open('out_file.pdb', 'r') as f:
                        st.download_button(label="Download file as PDB", data=f, file_name='out_file.pdb')

                
        elif file_name.endswith('pdb'):
            smiles_out = st.checkbox('SMILES', key=45)
            mol_out = st.checkbox('MOL', key=46)
            sdf_out = st.checkbox('SDF', key=47)
            if st.button("Convert"):
                with open('input_file_{}.pdb'.format(int(time.time())), 'w') as f:
                    f.write(uploaded_file.getvalue().decode("utf-8"))
                if smiles_out:
                    os.system("obabel input_file_{}.pdb -O out_file.smi".format(int(time.time())))

                    with open('out_file.smi', 'r') as f:
                        st.download_button(label="Download file as SMILES", data=f, file_name='out_file.smi')
                
                if mol_out:
                    os.system("obabel input_file_{}.sdf -O out_file.mol".format(int(time.time())))

                    with open('out_file.mol', 'r') as f:
                        st.download_button(label="Download file as SDF", data=f, file_name='out_file.mol')

                if sdf_out:
                    os.system("obabel input_file_{}.sdf -O out_file.pdb".format(int(time.time())))

                    with open('out_file.pdb', 'r') as f:
                        st.download_button(label="Download file as PDB", data=f, file_name='out_file.pdb')
        
        st.header('Protein conversions')
        col1, col2 = st.columns(2)


        with col1:
            st.subheader("Select input format")
            pdb_in = st.checkbox('PDB', key=9)
            xtc_in = st.checkbox('XTC', key=10)           

        with col2:
            st.subheader("Choose output format")
            pdb_out = st.checkbox('PDB', key=11)
            xtc_out = st.checkbox('XTC', key=12)
     
        if st.button("Convert", key=20):
            if pdb_in:
                with open('input_traj_{}.pdb'.format(int(time.time())), 'w') as f:
                    f.write(uploaded_file)
                if mol_out:
                    os.system("mdconvert input_traj_{}.pdb -o out_traj.xtc".format(int(time.time())))

                    with open('out_traj.xtc', 'r') as f:
                        st.download_button(label="Download file as xtc", data=f, file_name='out_traj.xtc')

            if xtc_in:
                with open('input_traj_{}.xtc'.format(int(time.time())), 'w') as f:
                    f.write(uploaded_file)
                if mol_out:
                    os.system("mdconvert input_traj_{}.xtc -o out_traj.pdb".format(int(time.time())))

                    with open('out_traj.pdb', 'r') as f:
                        st.download_button(label="Download file as pdb", data=f, file_name='out_traj.pdb')
            



if __name__=='__main__':
    main()
