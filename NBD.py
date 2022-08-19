import streamlit as st

def intro():
    import streamlit as st

    st.write("# Welcome to NBD Suite alpha! 👋")
    st.sidebar.success("Select a service above.")

    st.markdown(
        """
        **👈 Select a service from the dropdown on the left** to see some examples
        of what NBD Suite can do!

        ### Want to learn more?

        - Check out documentation [streamlit.io](https://nostrumbiodiscovery.github.io/nbd_central_docs/index.html)
        

        ### Got any feedback?
        Your comments and suggestions are really key for improving any of the products at any stage, so please, go ahead!
        - Push a ticket to our IT Team [here](https://forms.monday.com/forms/3d8b259d3d6269d9af3ebc976d51e9f4?r=use1)
        - Leave a rating and a suggestion [here](https://forms.gle/xkudVjsW2eMUAJxo6)
    """
    )

intro()