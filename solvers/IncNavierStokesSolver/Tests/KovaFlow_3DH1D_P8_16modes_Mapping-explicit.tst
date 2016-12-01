<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=8, 16 Fourier modes, using explicit mapping</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P8_16modes_Mapping-explicit.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P8_16modes_Mapping-explicit.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-9">1.46561e-05</value>
            <value variable="v" tolerance="1e-9">1.25542e-08</value>
            <value variable="w" tolerance="1e-9">1.08398e-07</value>
	    <value variable="p" tolerance="1e-9">1.08737e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9">1.64559e-05</value>
            <value variable="v" tolerance="1e-9">2.01086e-07</value>
            <value variable="w" tolerance="1e-9">3.48699e-07</value>
	    <value variable="p" tolerance="1e-9">8.01151e-06</value>
        </metric>
    </metrics>
</test>
