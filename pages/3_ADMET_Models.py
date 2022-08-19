import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors
from scopy.ScoDruglikeness import molproperty
import joblib
import streamlit.components.v1 as components


def compute_descriptors(input_file):

    in_df = pd.read_csv(input_file)
    desc_df = pd.DataFrame(columns=['canon_smiles', 'inchikey', 'mw', 'hba', 'hbd', 'rtb', 'heteroatoms',
    'aromatic_rings', 'heavy_atoms', 'psa', 'logp', 'logd', 'logsw', 'rgb', 'C',
    'n_stereocenters', 'Fsp3'])

    for i in in_df['smiles']:
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

def prediction(descriptor_file, model):

    if model != 'Rat acute toxicity':

        if model == 'Blood-brain barrier permeability':
            pickle_in = open('BBpred/blood_brain_barrier.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'CaCo2 permeability':
            pickle_in = open('CaCo/caco.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'Human intestinal absorption':
            pickle_in = open('hia/hia.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'cyp450 1a2 inhibition':
            pickle_in = open('cyp1a2/cyp1a2.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'cyp450 2c19 inhibition':
            pickle_in = open('cyp2c19/cyp2c19.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'cyp450 2c9 inhibition':
            pickle_in = open('cyp2c9/cyp2c9.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'cyp450 2d6 inhibition':
            pickle_in = open('cyp2d6/cyp2d6.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'cyp450 2d6 substrate':
            pickle_in = open('cyp2d6s/cyp2d6s.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'cyp450 3a4 inhibition':
            pickle_in = open('cyp3a4/cyp3a4.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'cyp450 3a4 substrate':
            pickle_in = open('cyp3a4s/cyp3a4s.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'cyp450 inhibition promiscuity':
            pickle_in = open('cypip/cypip.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'herg inhibition pred I':
            pickle_in = open('hergi1/hergi1.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'herg inhibition pred II':
            pickle_in = open('hergi2/hergi2.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'p glycoprotein inhibition I':
            pickle_in = open('pglyi1/pglyi1.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'p glycoprotein inhibition II':
            pickle_in = open('pglyi2/pglyi2.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'p glycoprotein substrate':
            pickle_in = open('pglys/pglys.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'Ames test':
            pickle_in = open('ames/ames.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'Biodegradation':
            pickle_in = open('biodeg/biodeg.pkl', 'rb')
            classifier = joblib.load(pickle_in)

        elif model == 'Carcinogenicity':
            pickle_in = open('carcino/carcino.pkl', 'rb')
            classifier = joblib.load(pickle_in)


        prediction = classifier.predict_proba(descriptor_file)

    else:
        pickle_in = open('ratac/ratac.pkl', 'rb')
        regressor = joblib.load(pickle_in)

        prediction = regressor.predict(descriptor_file)

    return prediction

def init_outdf(desc_df):
    out_df = pd.DataFrame(columns=['Compound'])
    out_df['Compound'] = desc_df['canon_smiles']

    return out_df

def build_prediction_table(out_df, desc_df, prediction, model):

    if model != 'Rat acute toxicity':
        probability_list = [round(x[0], 3) for x in prediction]
        out_df['{} probability'.format(model)] = probability_list
    else:
        metric_list = list(prediction)
        out_df['{}'.format(model)] = metric_list

    return out_df

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
    <h1 style ="color:white;text-align:center;">NBD ADMET Prediction Models</h1>
    </div>
    <div style ="background-color:white;padding:40px">
    </div>
    """


    sidebar_html = '''
    <h6 style ="color:Grey;text-align:center;">Supported models:</h6>
    <ul style="list-style-type:circle;">
    <li>Blood-brain barrier permeability</li>
    <li>CaCo2 permeability</li>
    <li>Human intestinal absorption</li>
    <li>Carcinogenicity</li>
    <li>cyp450 1a2 inhibition</li>
    <li>cyp450 2c19 inhibition</li>
    <li>cyp450 2c9 inhibition</li>
    <li>cyp450 2d6 inhibition</li>
    <li>cyp450 2d6 substrate</li>
    <li>cyp450 3a4 inhibition</li>
    <li>cyp450 3a4 substrate</li>
    <li>cyp450 inhibition promiscuity</li>
    <li>herg inhibition pred I</li>
    <li>herg inhibition pred II</li>
    <li>p glycoprotein inhibition I</li>
    <li>p glycoprotein inhibition II</li>
    <li>p glycoprotein substrate</li>
    <li>Rat acute toxicity</li>
    <li>Ames Test</li>
    <li>Biodegradation</li>
    </ul>
    <p><span class="small">v0.1. - 032022</span></p>
    '''

    st.markdown(html_temp, unsafe_allow_html = True)
    add_selectbox = st.sidebar.markdown(sidebar_html, unsafe_allow_html = True)

    uploaded_file = st.file_uploader("Choose a file")
    result =""

    options = st.multiselect(
     'Properties to predict - You may choose more than one model',
     ['Blood-brain barrier permeability', 'CaCo2 permeability', 'Human intestinal absorption',
      'Carcinogenicity', 'cyp450 1a2 inhibition','cyp450 2c19 inhibition',
      'cyp450 2c9 inhibition',  'cyp450 2d6 inhibition',
      'cyp450 2d6 substrate', 'cyp450 3a4 inhibition', 'cyp450 3a4 substrate',
      'cyp450 inhibition promiscuity', 'herg inhibition pred I',
      'herg inhibition pred II', 'p glycoprotein inhibition I',
      'p glycoprotein inhibition II', 'p glycoprotein substrate', 'Rat acute toxicity',
      'Ames Test', 'Biodegradation'])

    if uploaded_file is not None:

        if st.button("Predict"):
            with st.spinner('Computing descriptors...'):
                descr_file = compute_descriptors(uploaded_file)
            st.success('Done! Now fitting results...')
            out_df = init_outdf(descr_file)

            if len(options) > 0:

                for model in options:
                    
                    result = prediction(descr_file, model)
                    pred_df = build_prediction_table(out_df, descr_file, result, model)

                st.dataframe(pred_df)

                pred_df = convert_df(pred_df)
                descr_file  = convert_df(descr_file)

                st.download_button('Download prediction as CSV', data=pred_df, file_name='admet_pred.csv')

            else:
                 st.error('ERROR! No model selected. Please, select at least one model.')

    



if __name__=='__main__':
    main()
