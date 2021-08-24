<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 3D CG, hexes, mixed BCs, iterative SC, contiguous storage strategy</description>
    <executable>Helmholtz3D</executable>
    <parameters>-I GlobalSysSoln=IterativeStaticCond -I LocalMatrixStorageStrategy=Contiguous Helmholtz3D_Hex_AllBCs_P6.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Hex_AllBCs_P6.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-8">0.000416575</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-8">0.000871589</value>
        </metric>
    </metrics>
</test>


