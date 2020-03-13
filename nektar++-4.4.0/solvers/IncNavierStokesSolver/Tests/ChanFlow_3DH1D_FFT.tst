<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Laminar Channel Flow 3D homogeneous 1D, P=3, 20 Fourier modes (FFT)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_3DH1D_FFT.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_3DH1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">3.61998e-16</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">2.815e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">1.88738e-15</value>
            <value variable="v" tolerance="1e-6">2.78215e-16</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">1.4988e-13</value>
        </metric>
    </metrics>
</test>


