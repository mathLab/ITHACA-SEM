<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=10</description>
    <executable>APESolver</executable>
    <parameters>APE_2DPulseAdv_Galerkin_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_2DPulseAdv_Galerkin_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-12">21.3342</value>
            <value variable="u" tolerance="1e-12">0.0101474</value>
            <value variable="v" tolerance="1e-12">0.010618</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-12">93.234</value>
            <value variable="u" tolerance="1e-12">0.0384036</value>
            <value variable="v" tolerance="1e-12">0.0383323</value>
        </metric>
    </metrics>
</test>
