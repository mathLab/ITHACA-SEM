<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=10</description>
    <executable>APESolver</executable>
    <parameters>APE_2DPulseAdv_FRDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_2DPulseAdv_FRDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-12">21.3359</value>
            <value variable="u" tolerance="1e-12">0.0101433</value>
            <value variable="v" tolerance="1e-12">0.0106118</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-12">91.4563</value>
            <value variable="u" tolerance="1e-12">0.0382887</value>
            <value variable="v" tolerance="1e-12">0.0381487</value>
        </metric>
    </metrics>
</test>
