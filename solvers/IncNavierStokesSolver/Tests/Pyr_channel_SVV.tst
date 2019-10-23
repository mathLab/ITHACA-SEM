<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, Pyramidic elements, using SVV</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Pyr_channel_SVV.xml</parameters>
    <files>
        <file description="Session File">Pyr_channel_SVV.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10">2.67985e-10</value>
            <value variable="v" tolerance="1e-10">2.11332e-10</value>
            <value variable="w" tolerance="1e-9">8.34232e-10</value>
            <value variable="p" tolerance="1e-8">1.51191e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9">1.09756e-09</value>
            <value variable="v" tolerance="1e-9">9.06097e-10</value>
            <value variable="w" tolerance="1e-8">5.09333e-09</value>
            <value variable="p" tolerance="1e-6">3.10755e-07</value>
        </metric>
    </metrics>
</test>
