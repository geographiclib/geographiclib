using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Data;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using NETGeographicLib;

namespace Projections
{
    public partial class RhumbPanel : UserControl
    {
        Rhumb m_rhumb;
        public RhumbPanel()
        {
            InitializeComponent();
            m_rhumb = Rhumb.WGS84();
            m_Lat1TextBox.Text = "32";
            m_lon1TextBox.Text = "-86";
            m_lat2TextBox.Text = "33";
            m_lon2TextBox.Text = "-87";
            OnIndirect(this, null);
        }

        private void OnDirect(object sender, EventArgs e)
        {
            try
            {
                double lat1 = Double.Parse( m_Lat1TextBox.Text );
                double lon1 = Double.Parse( m_lon1TextBox.Text );
                double s12 = Double.Parse( m_s12TextBox.Text );
                double azimuth = Double.Parse( m_azimuthTextBox.Text );
                double lat2, lon2;
                m_rhumb.Direct(lat1, lon1, azimuth, s12, out lat2, out lon2);
                m_lat2TextBox.Text = lat2.ToString();
                m_lon2TextBox.Text = lon2.ToString();
                GeneratePoints(lat1, lon1, azimuth, 0.25 * s12);
            }
            catch (Exception xcpt )
            {
                MessageBox.Show( xcpt.Message, "Invalid input", MessageBoxButtons.OK, MessageBoxIcon.Error );
            }
        }

        private void OnIndirect(object sender, EventArgs e)
        {
            try
            {
                double lat1 = Double.Parse(m_Lat1TextBox.Text);
                double lon1 = Double.Parse(m_lon1TextBox.Text);
                double lat2 = Double.Parse(m_lat2TextBox.Text);
                double lon2 = Double.Parse(m_lon2TextBox.Text);
                double s12, azimuth;
                m_rhumb.Inverse(lat1, lon1, lat2, lon2, out s12, out azimuth);
                m_s12TextBox.Text = s12.ToString();
                m_azimuthTextBox.Text = azimuth.ToString();
                GeneratePoints(lat1, lon1, azimuth, 0.25 * s12);
            }
            catch (Exception xcpt)
            {
                MessageBox.Show(xcpt.Message, "Invalid input", MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
        }

        private void GeneratePoints( double lat1, double lon1, double azimuth, double space )
        {
            RhumbLine line = m_rhumb.Line(lat1, lon1, azimuth);
            m_pointsView.Items.Clear();
            for (int i = 0; i < 5; i++)
            {
                double lat2, lon2;
                line.Position(i * space, out lat2, out lon2);
                string[] items = new string[2] { lat2.ToString(), lon2.ToString() };
                ListViewItem item = new ListViewItem(items);
                m_pointsView.Items.Add(item);
            }
        }
    }
}
