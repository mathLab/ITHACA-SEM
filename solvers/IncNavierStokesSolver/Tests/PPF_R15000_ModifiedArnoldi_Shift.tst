<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linear stability with coupled solver (LM with Modified Arnoldiand Complex Shift): ChannelMax Ev = (2.4868e-03,1.5835e-01i) </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters> -I Driver=ModifiedArnoldi PPF_R15000_3D.xml</parameters>
    <files>
        <file description="Session File">PPF_R15000_3D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="2">
            <value variable="u" tolerance="0.0001">2.31851</value>
            <value variable="v" tolerance="0.0001">0.131819</value>
            <value variable="w" tolerance="0.0001">2.19656</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="u" tolerance="0.0001">0.79388</value>
            <value variable="v" tolerance="0.0001">0.0189237</value>
            <value variable="w" tolerance="0.0001">0.710008</value>
        </metric>
    </metrics>
</test>


