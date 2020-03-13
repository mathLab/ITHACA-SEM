<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=10</description>
    <executable>APESolver</executable>
    <parameters>APE_2DVariableC_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_2DVariableC_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-4">145.343</value>
            <value variable="u" tolerance="1e-7">0.329227</value>
            <value variable="v" tolerance="1e-7">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-4">119.917</value>
            <value variable="u" tolerance="1e-7">0.260865</value>
            <value variable="v" tolerance="1e-7">0</value>
        </metric>
    </metrics>
</test>
