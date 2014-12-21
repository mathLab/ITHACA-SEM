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
            <value variable="u">55.4979</value>
            <value variable="v">3.73647</value>
            <value variable="w">0.0065555</value>
	    	<value variable="p">166.89</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u">1.44445</value>
            <value variable="v">0.671072</value>
            <value variable="w">0.00443441</value>
	    	<value variable="p">5.74388</value>
        </metric>
    </metrics>
</test>
