<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D discontinuous material properties, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>VariableMatPropTest.xml</parameters>
    <files>
        <file description="Session File">VariableMatPropTest.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-12">5.47723</value>
            <value variable="u" tolerance="1e-12">5.47899</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-12">1.00008</value>
            <value variable="u" tolerance="1e-12">1.04038</value>
        </metric>
    </metrics>
</test>
