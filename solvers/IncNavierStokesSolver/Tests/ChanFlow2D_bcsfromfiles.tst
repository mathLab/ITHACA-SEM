<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=5 Boundary Conditions from files</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow2D_bcsfromfiles.xml</parameters>
    <files>
        <file description="Session File">ChanFlow2D_bcsfromfiles.xml</file>
        <file description="Session File">ChanFlow2D_bcsfromfiles_u_1.bc</file>
	<file description="Session File">ChanFlow2D_bcsfromfiles_u_3.bc</file>
	<file description="Session File">ChanFlow2D_bcsfromfiles_v_3.bc</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">1.85224e-13</value>
            <value variable="v" tolerance="1e-6">2.65241e-13</value>
	    <value variable="p" tolerance="1e-6">4.8272e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">8.31335e-13</value>
            <value variable="v" tolerance="1e-6">2.54448e-13</value>
	    <value variable="p" tolerance="1e-6">1.03323e-10</value>
        </metric>
    </metrics>
</test>


