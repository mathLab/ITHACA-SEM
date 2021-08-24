<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 2D CG with P=7, all BCs, iterative sc, par(2), check parallel dir bnd</description>
    <executable>Helmholtz2D</executable>
    <parameters>--use-scotch -I GlobalSysSoln=IterativeStaticCond Helmholtz2D_P7_AllBCs.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">Helmholtz2D_P7_AllBCs.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-7">0.00888406</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">0.0102587</value>
        </metric>
    </metrics>
</test>


