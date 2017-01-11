<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=3 using Weak Pressure form of VCS</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_m3_VCSWeakPress.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_m3_VCSWeakPress.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">5.92258e-16</value>
            <value variable="v" tolerance="1e-12">1.72565e-16</value>
	    <value variable="p" tolerance="1e-12">5.80063e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">4.30211e-15</value>
            <value variable="v" tolerance="1e-12">3.56426e-16</value>
	    <value variable="p" tolerance="1e-12">1.20737e-14</value>
        </metric>
    </metrics>
</test>


