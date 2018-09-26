<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=5, 20 Fourier modes - Skew-Symmetric advection(MVM)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P5_20modes_SKS_MVM.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P5_20modes_SKS_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.91624e-06</value>
            <value variable="v" tolerance="1e-12">1.42505e-06</value>
            <value variable="w" tolerance="1e-12">8.73217e-07</value>
	    <value variable="p" tolerance="1e-12">2.0521e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">7.40058e-06</value>
            <value variable="v" tolerance="1e-12">2.2922e-06</value>
            <value variable="w" tolerance="1e-12">1.94939e-06</value>
	    <value variable="p" tolerance="1e-12">6.07939e-05</value>
        </metric>
    </metrics>
</test>
