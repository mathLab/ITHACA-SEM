<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linear stability with coupled solver (LM with Modified Arnoldiand Complex Shift): ChannelMax Ev = (2.4868e-03,1.5835e-01i) </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters> -I Driver=ModifiedArnoldi PPF_R15000_3D.xml</parameters>
    <files>
        <file description="Session File">PPF_R15000_3D.xml</file>
    </files>
    <metrics>
        <metric type="Eigenvalue" id="0">
            <value index="0" tolerance="0.001">2.8057e-01,-2.4005e+01</value>
            <value index="1" tolerance="0.001">2.8057e-01,2.4005e+01</value>
        </metric>
    </metrics>
</test>


