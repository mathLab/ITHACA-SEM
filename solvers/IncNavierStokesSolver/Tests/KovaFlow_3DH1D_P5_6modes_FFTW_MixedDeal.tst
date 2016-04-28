<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=5, 6 Fourier modes (FFTW), specHP+homog dealiasing</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P5_6modes_FFTW_MixedDeal.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P5_6modes_FFTW_MixedDeal.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10">1.55951e-06</value>
            <value variable="v" tolerance="1e-10">9.42128e-07</value>
            <value variable="w" tolerance="1e-10">9.13083e-07</value>
	    <value variable="p" tolerance="1e-10">2.03518e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10">2.68238e-06</value>
            <value variable="v" tolerance="1e-10">2.08413e-06</value>
            <value variable="w" tolerance="1e-10">1.59952e-06</value>
	    <value variable="p" tolerance="1e-10">6.95616e-05</value>
        </metric>
    </metrics>
</test>
