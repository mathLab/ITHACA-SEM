<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Womersley B.C. Test</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Womersley_PipeFlow.xml</parameters>
    <files>
        <file description="Session File">Womersley_PipeFlow.xml</file>
        <file description="Session File">WomParams.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">0.00192525</value>
            <value variable="v" tolerance="1e-8">0.00171326</value>
	        <value variable="w" tolerance="1e-8">0.0338466</value>
	        <value variable="p" tolerance="1e-8">0.0740287</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">0.0102787</value>
            <value variable="v" tolerance="1e-8">0.00891606</value>
            <value variable="w" tolerance="1e-8">0.131224</value>
            <value variable="p" tolerance="1e-8">0.142029</value>
        </metric>
    </metrics>
</test>

