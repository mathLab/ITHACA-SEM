<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D tapering vessel, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>VariableAreaTest.xml</parameters>
    <files>
        <file description="Session File">VariableAreaTest.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-5">10.7976</value>
            <value variable="u" tolerance="1e-5">19.9526</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-5">0.988339</value>
            <value variable="u" tolerance="3e-4">1.99782</value>
        </metric>
    </metrics>
</test>


