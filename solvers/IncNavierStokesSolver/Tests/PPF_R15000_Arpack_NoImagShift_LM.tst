<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linear stability with coupled solver (LM with Arpack and Real Shift): ChannelMax Ev = (0.00248682 -0.158348i) </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters> -P nvec=4 -P kdim=32 -P imagShift=0.0 -I ArpackProblemType=LargestMag PPF_R15000_3D.xml</parameters>
    <files>
        <file description="Session File">PPF_R15000_3D.xml</file>
    </files>
    <metrics>
        <metric type="Eigenvalue" id="0">
            <value tolerance="0.001">-0.000201531,0</value>
        </metric>
    </metrics>
</test>


