<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linear stability with coupled solver (LI with Arpack and Real Shift): ChannelMax Ev = (0.00248682 -0.158348i) </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters> -P nvec=20 -P kdim=256 -P imagShift=0.0 -I ArpackProblemType=LargestImag PPF_R15000_3D.xml</parameters>
    <files>
        <file description="Session File">PPF_R15000_3D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="2">
            <value variable="u" tolerance="0.0001">0.481905</value>
            <value variable="v" tolerance="0.0001">0.231467</value>
            <value variable="w" tolerance="0.0001">0.628925</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="u" tolerance="0.0001">0.0901151</value>
            <value variable="v" tolerance="0.0001">0.0507247</value>
            <value variable="w" tolerance="0.0001">0.117699</value>
        </metric>
    </metrics>
</test>


