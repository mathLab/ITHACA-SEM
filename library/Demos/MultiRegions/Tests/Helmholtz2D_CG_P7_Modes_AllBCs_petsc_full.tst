<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz 2D CG with P=7, PETSc full</description>
    <executable>Helmholtz2D</executable>
    <parameters>-I GlobalSysSoln=PETScFull Helmholtz2D_P7_AllBCs.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_P7_AllBCs.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.047035</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">0.0959695</value>
        </metric>
    </metrics>
</test>
