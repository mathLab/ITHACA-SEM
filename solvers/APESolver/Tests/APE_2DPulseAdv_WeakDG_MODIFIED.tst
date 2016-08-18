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
            <value variable="p" tolerance="1e-4">21.3387</value>
            <value variable="u" tolerance="1e-7">0.0101445</value>
            <value variable="v" tolerance="1e-7">0.0105278</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-4">91.6371</value>
            <value variable="u" tolerance="1e-7">0.0383799</value>
            <value variable="v" tolerance="1e-7">0.0381269</value>
        </metric>
    </metrics>
</test>
