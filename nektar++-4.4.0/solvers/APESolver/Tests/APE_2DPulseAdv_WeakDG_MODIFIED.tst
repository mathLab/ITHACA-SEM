<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=10</description>
    <executable>APESolver</executable>
    <parameters>APE_2DPulseAdv_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_2DPulseAdv_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-4">6.77218</value>
            <value variable="u" tolerance="1e-7">0.00160109</value>
            <value variable="v" tolerance="1e-7">0.00163104</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-4">30.1206</value>
            <value variable="u" tolerance="1e-7">0.00609606</value>
            <value variable="v" tolerance="1e-7">0.00607659</value>
        </metric>
    </metrics>
</test>
