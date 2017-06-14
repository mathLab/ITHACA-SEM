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
            <value variable="u" tolerance="1e-10">2.6682e-10</value>
            <value variable="v" tolerance="1e-10">2.67015e-10</value>
            <value variable="w" tolerance="1e-9">1.15296e-09</value>
            <value variable="p" tolerance="1e-8">2.8591e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9">1.838e-09</value>
            <value variable="v" tolerance="1e-9">1.90318e-09</value>
            <value variable="w" tolerance="1e-8">1.42728e-08</value>
            <value variable="p" tolerance="1e-6">1.0353e-06</value>
        </metric>
    </metrics>

</test>
