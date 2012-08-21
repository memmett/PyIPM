# 
# Dictionary of mesh sizes used in the numerical experiments.
# 

mesh_sizes = {

    ('Exact', 'MidPoint'): [ 100, 200, 400, 800, 1600 ],
    ('Exact', 'Zuidema'):  [ 50, 100, 200 ],

    ('Exact', 'GENClarkQP'):     [ 100, 200 ],
    ('Exact', 'GENClarkGL(3)'):  [ 100, 200, 400 ],
    ('Exact', 'GENClarkGL(5)'):  [ 100, 200, 400 ],
    ('Exact', 'GENClarkGL(9)'):  [ 100, 200, 400 ],
    ('Exact', 'GENClarkCC(3)'):  [ 100, 200, 400 ],
    ('Exact', 'GENClarkCC(5)'):  [ 100, 200, 400 ],
    ('Exact', 'GENClarkCC(9)'):  [ 100, 200, 400 ],

    ('Exact', 'INTClarkQP'):     [ 100, 200 ],
    ('Exact', 'INTClarkGL(3)'):  [ 100, 200, 400 ],
    ('Exact', 'INTClarkGL(5)'):  [ 100, 200, 400 ],
    ('Exact', 'INTClarkGL(9)'):  [ 100, 200, 400 ],
    ('Exact', 'INTClarkCC(3)'):  [ 100, 200, 400 ],
    ('Exact', 'INTClarkCC(5)'):  [ 100, 200, 400 ],
    ('Exact', 'INTClarkCC(9)'):  [ 100, 200, 400 ],

    ('Exact', 'GaussLegendre'):    [ 20, 40, 80, 160, 320, 640 ],
    ('Exact', 'GaussLegendre(3)'): [ 50, 100, 200, 400 ],
    ('Exact', 'GaussLegendre(5)'): [ 50, 100, 200, 400 ],
    ('Exact', 'GaussLegendre(9)'): [ 50, 100, 200, 400 ],

    ('Exact', 'ClenshawCurtis'):    [ 20, 40, 80, 160, 320, 640 ],
    ('Exact', 'ClenshawCurtis(3)'): [ 50, 100, 200, 400 ],
    ('Exact', 'ClenshawCurtis(5)'): [ 50, 100, 200, 400 ],
    ('Exact', 'ClenshawCurtis(9)'): [ 50, 100, 200, 400 ],

    ('EED', 'MidPoint'): [ 100, 200, 400, 800, 1600 ],
    ('EED', 'Zuidema'):  [ 50, 100, 200, 400 ],

    ('EED', 'GENClarkQP'):     [ 100, 200 ],
    ('EED', 'GENClarkGL(3)'):  [ 100, 200, 400 ],
    ('EED', 'GENClarkGL(5)'):  [ 100, 200, 400 ],
    ('EED', 'GENClarkGL(9)'):  [ 100, 200, 400 ],
    ('EED', 'GENClarkCC(3)'):  [ 100, 200, 400 ],
    ('EED', 'GENClarkCC(5)'):  [ 100, 200, 400 ],
    ('EED', 'GENClarkCC(9)'):  [ 100, 200, 400 ],

    ('EED', 'INTClarkQP'):     [ 100, 200 ],
    ('EED', 'INTClarkGL(3)'):  [ 100, 200, 400 ],
    ('EED', 'INTClarkGL(5)'):  [ 100, 200, 400 ],
    ('EED', 'INTClarkGL(9)'):  [ 100, 200, 400 ],
    ('EED', 'INTClarkCC(3)'):  [ 100, 200, 400 ],
    ('EED', 'INTClarkCC(5)'):  [ 100, 200, 400 ],
    ('EED', 'INTClarkCC(9)'):  [ 100, 200, 400 ],

    ('EED', 'GaussLegendre'):    [ 20, 40, 80, 160, 320, 640 ],
    ('EED', 'GaussLegendre(3)'): [ 50, 100, 200, 400 ],
    ('EED', 'GaussLegendre(5)'): [ 50, 100, 200, 400 ],
    ('EED', 'GaussLegendre(9)'): [ 50, 100, 200, 400 ],

    ('EED', 'AdjGaussLegendre(3)'): [ 50, 100, 200, 400 ],
    ('EED', 'AdjGaussLegendre(5)'): [ 50, 100, 200, 400 ],
    ('EED', 'AdjGaussLegendre(9)'): [ 50, 100, 200, 400 ],

    ('EED', 'ClenshawCurtis'):    [ 20, 40, 80, 160, 320, 640 ],
    ('EED', 'ClenshawCurtis(3)'): [ 50, 100, 200, 400 ],
    ('EED', 'ClenshawCurtis(5)'): [ 50, 100, 200, 400 ],
    ('EED', 'ClenshawCurtis(9)'): [ 50, 100, 200, 400 ],

    ('Zuidema', 'MidPoint'): [ 200, 400, 800, 1600 ],
    ('Zuidema', 'Zuidema'):  [ 200, 400, 800 ],

    ('Zuidema', 'GENClarkQP'):     [ 200, 400 ],
    ('Zuidema', 'GENClarkGL(3)'):  [ 200, 400, 800 ],
    ('Zuidema', 'GENClarkGL(5)'):  [ 200, 400, 800 ],
    ('Zuidema', 'GENClarkGL(9)'):  [ 200, 400, 800 ],
    ('Zuidema', 'GENClarkCC(3)'):  [ 200, 400, 800 ],
    ('Zuidema', 'GENClarkCC(5)'):  [ 200, 400, 800 ],
    ('Zuidema', 'GENClarkCC(9)'):  [ 200, 400, 800 ],

    ('Zuidema', 'INTClarkQP'):     [ 200, 400 ],
    ('Zuidema', 'INTClarkGL(3)'):  [ 200, 400, 800 ],
    ('Zuidema', 'INTClarkGL(5)'):  [ 200, 400, 800 ],
    ('Zuidema', 'INTClarkGL(9)'):  [ 200, 400, 800 ],

    ('Zuidema', 'GaussLegendre'):    [ 80, 160, 320, 640 ],
    ('Zuidema', 'GaussLegendre(3)'): [ 200, 400, 800 ],
    ('Zuidema', 'GaussLegendre(5)'): [ 200, 400, 800 ],
    ('Zuidema', 'GaussLegendre(9)'): [ 200, 400, 800 ],

    ('Zuidema', 'ClenshawCurtis'):    [ 80, 160, 320, 640 ],
    ('Zuidema', 'ClenshawCurtis(3)'): [ 200, 400, 800 ],
    ('Zuidema', 'ClenshawCurtis(5)'): [ 200, 400, 800 ],
    ('Zuidema', 'ClenshawCurtis(9)'): [ 200, 400, 800 ],

    }
