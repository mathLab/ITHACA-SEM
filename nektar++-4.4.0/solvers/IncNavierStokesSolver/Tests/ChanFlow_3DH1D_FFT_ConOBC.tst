<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Laminar Channel Flow 3D homogeneous 1D, P=3, 20 Fourier modes (FFT), Convective outflow</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_3DH1D_FFT_ConOBC.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_3DH1D_FFT_ConOBC.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-5">0</value>
            <value variable="v" tolerance="1e-5">0</value>
            <value variable="w" tolerance="1e-5">0</value>
            <value variable="p" tolerance="5e-5">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5">0</value>
            <value variable="v" tolerance="1e-5">0</value>
            <value variable="w" tolerance="1e-5">0</value>
            <value variable="p" tolerance="5e-5">0</value>
        </metric>
    </metrics>
</test>


