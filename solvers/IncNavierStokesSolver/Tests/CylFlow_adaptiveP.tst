<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D cylinder flow simulation using adaptive polynomial order</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>CylFlow_adaptiveP.xml</parameters>
    <files>
        <file description="Session File">CylFlow_adaptiveP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">21.4162</value>
            <value variable="v" tolerance="1e-12">0.911291</value>
            <value variable="p" tolerance="1e-12">0.875288</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.76892</value>
            <value variable="v" tolerance="1e-12">0.837794</value>
            <value variable="p" tolerance="1e-12">1.29043</value>
        </metric>
    </metrics>
</test>

