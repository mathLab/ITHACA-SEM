<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=5, 20 Fourier modes (MVM), SVV Homogeneous 1D</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P5_20modes_MVM_SVVHomo1D.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P5_20modes_MVM_SVVHomo1D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.54682e-06</value>
            <value variable="v" tolerance="1e-12">9.50615e-07</value>
            <value variable="w" tolerance="1e-12">9.03459e-07</value>
	    <value variable="p" tolerance="1e-12">2.04388e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">3.14479e-06</value>
            <value variable="v" tolerance="1e-12">2.25387e-06</value>
            <value variable="w" tolerance="1e-12">1.80791e-06</value>
	    <value variable="p" tolerance="1e-12">6.87018e-05</value>
        </metric>
    </metrics>
</test>
