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
            <value variable="v" tolerance="1e-12">3.73645</value>
            <value variable="w" tolerance="1e-12">0.000112658</value>
            <value variable="p" tolerance="1e-12">166.726</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.44443</value>
            <value variable="v" tolerance="1e-12">0.670852</value>
            <value variable="w" tolerance="1e-12">6.8828e-06</value>
            <value variable="p" tolerance="1e-12">5.61629</value>
        </metric>
    </metrics>
</test>
