<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linear stability with coupled solver (LI with Arpack and Real Shift): ChannelMax Ev = (0.00248682 -0.158348i) </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters> -P nvec=20 -P kdim=256 -P imagShift=0.0 -I ArpackProblemType=LargestImag PPF_R15000_3D.xml</parameters>
    <files>
        <file description="Session File">PPF_R15000_3D.xml</file>
    </files>
    <metrics>
        <metric type="Eigenvalue" id="0">
            <value index="0" tolerance="0.001">-0.016749,0</value>
        </metric>
    </metrics>
</test>


