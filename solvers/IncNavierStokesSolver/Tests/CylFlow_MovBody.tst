<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D flexible cylinder flow simulation using "MovingBody" module</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>CylFlow_MovBody.xml</parameters>
    <files>
        <file description="Session File">CylFlow_MovBody.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">55.4979</value>
            <value variable="v" tolerance="1e-12">3.73647</value>
            <value variable="w" tolerance="1e-12">0.0065555</value>
            <value variable="p" tolerance="1e-12">166.89</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.44445</value>
            <value variable="v" tolerance="1e-12">0.671072</value>
            <value variable="w" tolerance="1e-12">0.00443441</value>
            <value variable="p" tolerance="1e-12">5.74388</value>
        </metric>
    </metrics>
</test>
