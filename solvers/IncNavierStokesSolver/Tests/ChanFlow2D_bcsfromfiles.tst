<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=5 Boundary Conditions from files</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow2D_bcsfromfiles.xml</parameters>
    <files>
        <file description="Session File">ChanFlow2D_bcsfromfiles.xml</file>
        <file description="Session File">ChanFlow2D_bcsfromfiles_u1_0.bc</file>
	<file description="Session File">ChanFlow2D_bcsfromfiles_u3_0.bc</file>
	<file description="Session File">ChanFlow2D_bcsfromfiles_v3_0.bc</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">2.03944e-13</value>
            <value variable="v" tolerance="1e-6">1.30205e-13</value>
	    <value variable="p" tolerance="1e-6">5.35429e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">7.21201e-13</value>
            <value variable="v" tolerance="1e-6">1.18131e-13</value>
	    <value variable="p" tolerance="1e-6">1.03743e-10</value>
        </metric>
    </metrics>
</test>


