<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=10</description>
    <executable>AcousticSolver</executable>
    <parameters>APE_2DVariableC_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_2DVariableC_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-4">47.2682</value>
            <value variable="u" tolerance="1e-7">0.242422</value>
            <value variable="v" tolerance="1e-7">0.0671381</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-4">58.295</value>
            <value variable="u" tolerance="1e-7">0.278645</value>
            <value variable="v" tolerance="1e-7">0.105514</value>
        </metric>
    </metrics>
</test>

