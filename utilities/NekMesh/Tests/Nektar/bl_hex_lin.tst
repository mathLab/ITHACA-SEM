<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Nektar linear hexahedron with boundary layer splitting</description>
    <executable>NekMesh</executable>
    <parameters>-m bl:layers=3 -m jac:list bl_hex_lin.rea bl_hex_lin-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">bl_hex_lin.rea</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
