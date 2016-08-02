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
            <value variable="u" tolerance="5e-6">1.03944e-05</value>
            <value variable="v" tolerance="2e-6">4.98402e-06</value>
	    <value variable="p" tolerance="1e-6">2.85015e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="4e-5">6.46273e-05</value>
            <value variable="v" tolerance="1e-5">1.86564e-05</value>
	    <value variable="p" tolerance="1e-5">4.56732e-05</value>
        </metric>
    </metrics>
</test>


