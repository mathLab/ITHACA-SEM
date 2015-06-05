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
            <value variable="p" tolerance="1e-11">126.796</value>
            <value variable="u" tolerance="1e-11">0.308052</value>
            <value variable="v" tolerance="1e-11">3e-12</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-11">102.142</value>
            <value variable="u" tolerance="1e-11">0.2588</value>
            <value variable="v" tolerance="1e-11">3e-12</value>
        </metric>
    </metrics>
</test>
