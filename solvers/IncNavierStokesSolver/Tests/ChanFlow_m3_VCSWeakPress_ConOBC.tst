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
            <value variable="u" tolerance="1e-6">2.36927e-06</value>
            <value variable="v" tolerance="1e-6">1.32509e-06</value>
	    <value variable="p" tolerance="1e-5">2.183e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5">1.88785e-05</value>
            <value variable="v" tolerance="3.0e-6">4.06239e-06</value>
	    <value variable="p" tolerance="1e-5">3.45384e-05</value>
        </metric>
    </metrics>
</test>


