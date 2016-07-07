<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=3 using Weak Pressure form of VCS and convective link outflow BC</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_m3_VCSWeakPress_ConOBC.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_m3_VCSWeakPress_ConOBC.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.07242e-06</value>
            <value variable="v" tolerance="1e-12">6.54816e-07</value>
	    <value variable="p" tolerance="1e-12">2.48295e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">8.06669e-06</value>
            <value variable="v" tolerance="1e-12">2.20307e-06</value>
	    <value variable="p" tolerance="1e-12">3.21805e-05</value>
        </metric>
    </metrics>
</test>


