<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=5, 20 Fourier modes (MVM)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P5_20modes_MVM.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P5_20modes_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.93348e-06</value>
            <value variable="v" tolerance="1e-12">1.41911e-06</value>
            <value variable="w" tolerance="1e-12">8.39121e-07</value>
	    <value variable="p" tolerance="1e-12">2.04066e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">7.34071e-06</value>
            <value variable="v" tolerance="1e-12">2.30025e-06</value>
            <value variable="w" tolerance="1e-12">1.80249e-06</value>
	    <value variable="p" tolerance="1e-12">6.10466e-05</value>
        </metric>
    </metrics>
</test>
