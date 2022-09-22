Tutorial 2: graph visualization
===============================

This tutorial will show you how to turn your Pasmopy Text model into a graph and generate both static and dynamic images from it.  
For demonstration purposes we will use a text representation of the `nfkb_pathway <https://github.com/biomass-dev/biomass/tree/master/biomass/models/nfkb_pathway>`_ model included in biomass. For a detailed description of the model, please refer to the following paper:  

* Oppelt, A. *et al*. Model-based identification of TNFα-induced IKKβ-mediated and IκBα-mediated regulation of NFκB signal transduction as a tool to quantify the impact of drug-induced liver injury compounds. *npj Syst. Biol. Appl.* **4**, 23 (2018). https://doi.org/10.1038/s41540-018-0058-z

Requirements
------------
* ``biomass>=0.9.0``
* ``pygraphviz>=1.9``
* ``pyvis>=0.2.1``
* ``graphviz>=2.42`` Installation instructions can be found `here <https://graphviz.org/download/>`_

Prepare a text file describing the model
----------------------------------------

.. code-block::
    :linenos:
    
    TNF synthesizes TNFR| kf=1 | TNF=1
    TNFR is degraded | kf=0.001
    TNFR + Ikk --> pIkk + TNFR | kf=0.0714 | Ikk=1
    pIkk is phosphorylated --> ppIkk | kf=0.0648
    ppIkk --> iIkk | kf=0.166
    iIkk --> Ikk | kf=0.0041
    pIkk + NfkIkb --> NfkpIkb + pIkk | kf=0.398 | NfkIkb=1
    pNfkIkb is phosphorylated --> pNfkpIkb | kf=1.3897
    pIkk + pNfkIkb --> pNfkpIkb + pIkk | kf=0.389
    pIkk + NfkIkb --> pNfkIkb + pIkk | kf=0.6438
    pIkk + NfkpIkb --> pNfkpIkb + pIkk | kf=0.2816
    NfkpIkb -->  Nfk + pIkb | kf=0.0811
    pNfkpIkb --> pNfk + pIkb | kf=1
    Nfk + Ikb --> NfkIkb | kf=2.839
    nNfk synthesizes mIkb | kf=0.0047
    mIkb is degraded | kf=0.0313
    mIkb synthesizes Ikb | kf=1
    pIkb is degraded | kf=0.6308
    Ikb translocates from cytoplasm to nucleus (1, 1) --> nIkb | kf=0.1226
    pNfk translocates from cytoplasm to nucleus (1, 1) --> pnNfk | kf=0.179585
    Nfk translocates from cytoplasm to nucleus (1, 1) --> nNfk | kf=0.01
    pnNfk --> nNfk | kf=1000
    nIkb binds nNfk --> nNfkIkb | kf=1000
    nNfkIkb translocates from nucleus to cytoplasm (1, 1) --> NfkIkb | kf=1000
    nNfk synthesizes RnaA20_1 | kf=1
    RnaA20_1 --> RnaA20 | kf=0.0311
    RnaA20 is degraded | kf=0.0089
    RnaA20 synthesizes A20 | kf=0.0006
    A20 is degraded | kf=0.0116

    @obs nuclear_IkBa: u[nIkb]
    @obs nuclear_NFkB: u[nNfk]

    @sim tspan: [0, 200]


.. note::
    For further details about setting up a biomass model from text please refer to the corresponding `biomass <https://biomass-core.readthedocs.io
    /en/latest/tutorial/cfos.html>`_ and `pasmopy <https://pasmopy.readthedocs.io/en/latest/model_development.html>`_ tutorials.

Import the model
----------------
.. code-block:: python

    from biomass import Text2Model
    
    model = Text2Model('name_of_your_txt_file.txt')
    
    model.convert()

Graph generation and static image
---------------------------------
The graph is constructed from the kinetic information gained during model construction. Connections always go from reactants/modifiers to products. There is no destinction made between modifiers and reactants, as well as activating and inhibiting modifiers. 

.. code-block:: python

    model.graph
    
The graph property contains an instance of the AGraph class implemented by pygraphviz. For available methods please refer to `their documentation <https://pygraphviz.github.io/documentation/stable/reference/agraph.html>`_. You can for example manually add/remove nodes or save the graph into a .dot file and import it into another 3rd party software.

A static image of the graph is drawn using

.. code-block:: python

    model.static_plot(save_dir='example_dir', file_name='nfkb_static.png')
    model.static_plot(save_dir='example_dir', file_name='nfkb_static_cust.png',
                      gviz_args='-Nshape=parallelogram -Nstyle=bold -Estyle=dashed')
    
The desired file format is inferred from the ending of file_name. Graphviz provides a variety of different engines that automatically generate a layout for the graph. By default the 'dot' engine is used, since it uses a hierarchical approach that is natural for biological data. Feel free to play around with the available engines, but be aware that biological networks can quickly become messy due to the prevalance of feedback interactions.  
Additionally graphviz provides a large variety of customization options, that have to be passed in the command line format. For a comprehensive list see the `graphviz manual <https://graphviz.org/pdf/dot.1.pdf>`_.  

.. figure:: ../_static/img/static_nfkb_graph.svg
    :align: center

Dynamic image
--------------
Thanks to the package `pyvis <https://github.com/WestHealth/pyvis>`_ we can also provide an interactive graph. The generation is just as simple as for the static image:  

.. code-block:: python

    model.dynamic_plot(save_dir='example_dir', file_name='nfkb_dynamic.html' show_controls=True, which_controls=['physics', 'layout'])
    
By default the plot will be immediately displayed in your browser. Set ``show`` to :obj:`False` if you don't want that. ``pyvis`` provides a variety of customization options as well. They can be directly accessed in the html file by setting ``show_controls`` to :obj:`True`. You can also specify which controls you want.

.. raw:: html
    :file: dynamic_nfkb_graph.html
