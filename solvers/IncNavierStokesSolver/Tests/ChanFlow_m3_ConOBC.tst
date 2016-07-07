<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=3 with convective link outflow BC</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_m3_ConOBC.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_m3_ConOBC.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">5.31595e-06</value>
            <value variable="v" tolerance="1e-12">3.47264e-06</value>
	    <value variable="p" tolerance="1e-12">2.8014e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">3.28775e-05</value>
            <value variable="v" tolerance="1e-12">1.41615e-05</value>
	    <value variable="p" tolerance="1e-12">3.62358e-05</value>
        </metric>
    </metrics>
</test>


