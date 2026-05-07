ARCS 1.5.1 

`setup_functions.py`

- removed ApplyDataToReactions as this class was now defunct 

`ReactionGibbsandEquilibrium`
- implemented reactit output which is more pythonic and faster 
- docstrings added 
- code readability now much better 
- tests added for this class
- Gibbs Free Energy of reaction is now per reactant atom 
    - this helps with unreasonably large K values 
    - the cost function in `GraphGenerator` has now removed the per reactant atom division as a result

`GraphGenerator` 
- removed trange and prange as it is not necessary anymore 
- removed multiprocessing functions as they were not necessary anymore 
- much faster
- code readability much better
- added tests for this class

`errors.py` 
-this is a module that will can calculate errors based on the supplied quantum calculated data - currently not used or implemented in arcs (work-in-progress)

`arcs-dash-app` 
- removed dash app to https://github.com/badw/arcs-dash.git

`Coupled Cluster Data` 
- now more compounds have been added to the backend data - the total list currently available are shown in `src/arcs/data/README.md` 

`pytests`
- tests for `analysis.py` and `traversal.py` have now been added






